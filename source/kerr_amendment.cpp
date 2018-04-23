//
//  kerr_amendment.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

const std::string KerrAmendment::exeption_msg = "Integral by $VAR has not trustable at $TIME, $RHO, $PHI, $Z";

/* Nonlinear medium */

KerrMedium::KerrMedium (double realative_mu, double realative_eps, double kerr, double conduct) 
: Homogeneous(realative_mu, realative_eps) 
{
	if (kerr >= 0) this->xi3 = kerr;
	else std::invalid_argument("Kerr coefficient must be positive.");
	
	if (conduct > 0) this->sigma = conduct;
	else std::invalid_argument("Kerr coefficient must be positive.");
}

double KerrMedium::conductivity (double ct, double z) const
{
	UNUSED(ct); UNUSED(z);
	return this->sigma;
}

double KerrMedium::relative_permittivity (double ct, double z, std::size_t term) const
{
	UNUSED(ct); UNUSED(z);
	switch(term) {
		// case 0: return 0;
		case 1: return Homogeneous::relative_permittivity(ct,z);
		case 3: return this->xi3;
		default: return 0;
	}
}

double KerrMedium::relative_permeability (double ct, double z, std::size_t term) const
{
	UNUSED(ct); UNUSED(z);
	switch(term) {
		// case 0: return 0;
		case 1: return Homogeneous::relative_permeability(ct,z);
		case 3: return 0;
		default: return 0;
	}
}

/* Nonlinear ammendmend */

double KerrAmendment::current_x (double vt, double rho, double phi, double z) const
{
	const double current_rho = this->current_rho(vt,rho,phi,z);
	const double current_phi = this->current_phi(vt,rho,phi,z);
	return current_rho * std::cos(phi) - current_phi * std::sin(phi);
}

double KerrAmendment::current_rho (double vt, double rho, double phi, double z) const
{
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double em_relation = std::sqrt(NonlinearMedium::MU0 * mu_r / NonlinearMedium::EPS0 * eps_r);
	double A3 = this->A0 * this->A0 * this->A0;

	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	double mult1 = i1_perp * (i2 - i1);
	double mult2 = i1 * (i2_perp - i1_perp);

	double coeff = A3 * em_relation * em_relation * em_relation / 8;
	double cos3sin0 = 3 * i1 * i1 * i1_perp * std::cos(phi) * std::cos(phi) * std::cos(phi);
	double cos1sin2 = (i2 - i1) * (mult1 + 2 * mult2) * std::cos(phi) * std::sin(phi) * std::sin(phi);
	return coeff * (cos3sin0 + cos1sin2);
}

double KerrAmendment::current_phi (double vt, double rho, double phi, double z) const
{
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double em_relation = std::sqrt(NonlinearMedium::MU0 * mu_r / NonlinearMedium::EPS0 * eps_r);
	double A3 = this->A0 * this->A0 * this->A0;

	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	double mult1 = i1_perp * (i2 - i1);
	double mult2 = i1 * (i2_perp - i1_perp);

	double coeff = A3 * em_relation * em_relation * em_relation / 8;
	double cos0sin3 = 3 * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp) * std::sin(phi) * std::sin(phi) * std::sin(phi);
	double cos2sin1 = i1 * (2 * mult1 + mult2) * std::cos(phi) * std::cos(phi) * std::sin(phi);
	return - coeff * (cos0sin3 + cos2sin1);
}

KerrAmendment::KerrAmendment (MissileField* field, KerrMedium* medium, UniformPlainDisk* source, Logger* logger_ptr)
: NonlinearField(field, medium) 
{
	this->global_logger = logger_ptr;
	this->linear_field = field;
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
}

double KerrAmendment::electric_x (double vt, double rho, double phi, double z) const
{
	double vt_z = vt - z;
	if (z == 0) return 0;
	if (vt_z < 1e-9) return 0;

	Logger* log_ptr = this->global_logger;

	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double em_relation = NonlinearMedium::MU0 * mu_r / NonlinearMedium::EPS0 * eps_r;
	double R = this->R, A0 = this->A0; 
	double coeff = kerr * A0 * A0 * A0 * em_relation * em_relation / 128;

	double integral;

	try {

		integral = SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, 2*R,
			[log_ptr, R, vt, rho, phi, z] (double z_perp) {

				try {

					return SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, 2*R,
						[log_ptr, R, vt, rho, phi, z, z_perp] (double rho_perp) {

							double max_vt = vt - z + z_perp; // grater then zero
							double max_nu = PERIODS_NU * std::abs(rho - rho_perp);
							double period_nu =  NODES_NU * (rho + rho_perp);
							std::size_t nodes_nu = max_nu * period_nu;
							if (!nodes_nu) nodes_nu = 1;

							double delta_sum;

							try  {

								delta_sum = SimpsonRunge(nodes_nu, MAX_ERROR).value(0, max_nu,
									[R, rho, phi, max_vt, rho_perp, z_perp] (double nu) {

										double sum = 0;
										sum += KerrAmendment::x_trans( -1, nu, rho, phi) *
											   KerrAmendment::N_sum(R, -1, nu, max_vt, rho_perp, z_perp);
										sum += KerrAmendment::x_trans( 1, nu, rho, phi) *
											   KerrAmendment::N_sum(R, 1, nu, max_vt, rho_perp, z_perp);
										sum += KerrAmendment::x_trans( -3, nu, rho, phi) *
											   KerrAmendment::N_sum(R, -3, nu, max_vt, rho_perp, z_perp);
										sum += KerrAmendment::x_trans( 3, nu, rho, phi) *
											   KerrAmendment::N_sum(R, 3, nu, max_vt, rho_perp, z_perp);
										return sum;
								} );	// <--------------------------------------------------------------------------------------- d nu1

							} catch (double not_trusted) {

								if (log_ptr) {
									std::string mesg = KerrAmendment::exeption_msg;
									mesg = std::regex_replace(mesg, std::regex("\\$VAR" ), "Z'");
									mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
									mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
									mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(phi));
									mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
									mesg += std::string(" (Z'= ") + std::to_string(z_perp);
									mesg += std::string(", RHO'= ") + std::to_string(rho_perp);
									mesg += std::string(")");
									log_ptr->warning(mesg);
								}

								delta_sum = not_trusted;
							}
							
							if (std::isnan(delta_sum)) delta_sum = 0;

							double step_sum = 0;
							try {

								step_sum = SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, max_vt,
									[log_ptr, R, vt, rho, phi, z, rho_perp, z_perp] (double vt_perp) {

									double delta_vt = vt - vt_perp;
									double delta_z = z - z_perp;
									double casual = std::sqrt(delta_vt*delta_vt - delta_z*delta_z);
									if (std::isnan(casual)) return 0.0;
									double period_nu = NODES_NU * std::abs(rho + rho_perp + casual);
									double max_nu =  PERIODS_NU * std::abs(rho - rho_perp + casual);
									std::size_t nodes_nu = max_nu * period_nu;
									if (!nodes_nu) nodes_nu = 1;

									double res;
									try {

										res = SimpsonRunge(nodes_nu, MAX_ERROR).value( 0, max_nu,
											[R, rho, phi, rho_perp, z_perp, vt_perp, casual, delta_vt] (double nu) {
												double bessel = jn(0,nu * casual) + jn(2,nu * casual); 
												double sum = 0;
												sum += KerrAmendment::x_trans( -1, nu, rho, phi) *
													   KerrAmendment::N_sum(R, -1, nu, vt_perp, rho_perp, z_perp);
												sum += KerrAmendment::x_trans( 1, nu, rho, phi) *
													   KerrAmendment::N_sum(R, 1, nu, vt_perp, rho_perp, z_perp);
												sum += KerrAmendment::x_trans( -3, nu, rho, phi) *
													   KerrAmendment::N_sum(R, -3, nu, vt_perp, rho_perp, z_perp);
												sum += KerrAmendment::x_trans( 3, nu, rho, phi) *
													   KerrAmendment::N_sum(R, 3, nu, vt_perp, rho_perp, z_perp);
												return sum * nu * nu * delta_vt * bessel / 2;

										} );	// <-------------------------------------------------------------------------------- d nu2

									} catch (double not_trusted) {

										if (log_ptr) {
											std::string mesg = KerrAmendment::exeption_msg;
											mesg = std::regex_replace(mesg, std::regex("\\$VAR" ), "Z'");
											mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
											mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
											mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(phi));
											mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
											mesg += std::string(" (Z'= ") + std::to_string(z_perp);
											mesg += std::string(", RHO'= ") + std::to_string(rho_perp);
											mesg += std::string(", VT'= ") + std::to_string(vt_perp);
											mesg += std::string(")");
											log_ptr->warning(mesg);
										}

										res = not_trusted;
									}

									if (std::isnan(res)) res = 0;
									return res;

								} );	// <---------------------------------------------------------------------------------------- d tau'

							} catch (double not_trusted) {

								if (log_ptr) {
									std::string mesg = KerrAmendment::exeption_msg;
									mesg = std::regex_replace(mesg, std::regex("\\$VAR" ), "Z'");
									mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
									mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
									mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(phi));
									mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
									mesg += std::string(" (Z'= ") + std::to_string(z_perp);
									mesg += std::string(", RHO'= ") + std::to_string(rho_perp);
									mesg += std::string(")");
									log_ptr->warning(mesg);
								}

								step_sum = not_trusted;
							}

						return delta_sum - step_sum;

					} ); //	<------------------------------------------------------------------------------------------------------- d rho'

				} catch (double not_trusted) {

					if (log_ptr) {
						std::string mesg = KerrAmendment::exeption_msg;
						mesg = std::regex_replace(mesg, std::regex("\\$VAR" ), "Z'");
						mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
						mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
						mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(phi));
						mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
						mesg += std::string(" (Z'= ") + std::to_string(z_perp);
						mesg += std::string(")");
						log_ptr->warning(mesg);
					}

					return not_trusted;
				}
		} );	// <---------------------------------------------------------------------------------------------------------------- d z'

	} catch (double not_trusted) {

		if (log_ptr) {
			std::string mesg = KerrAmendment::exeption_msg;
			mesg = std::regex_replace(mesg, std::regex("\\$VAR" ), "Z'");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(phi));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			log_ptr->warning(mesg);
		}

		integral = not_trusted;

		std::cout << "Ex (" << vt << ',' << rho << ',' << phi << ',' << z << ") => ";
		std::cout << - coeff * integral << std::endl;
	}

	return - coeff * integral;
}

double KerrAmendment::electric_rho (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::electric_phi (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::electric_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

double KerrAmendment::magnetic_rho (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

double KerrAmendment::x_trans (int m, double nu, double rho, double phi)
{
	double vect_rho = std::cos(phi) * std::cos(m * phi);
	double vect_phi = std::sin(phi) * std::sin(m * phi);
	vect_rho *= jn(m-1, nu * rho) + jn(m+1, nu * rho);
	vect_phi *= jn(m-1, nu * rho) - jn(m+1, nu * rho);
	return vect_rho + vect_phi;
}

double KerrAmendment::N_sum (double R, int m, double nu, double vt, double rho, double z)
{
	if ( vt   == 0 ) return 0.0;
	if ( z    == 0 ) return 0.0;
	if ( nu   == 0 ) return 0.0; // becouse math
	if ( rho  == 0 ) return 0.0; // becouse derivative
	if ( vt-z <= 0 ) return 0.0; // casuality princip
	double terms = 0;

	switch (m) {
		case -1: { 
			terms += 3 * KerrAmendment::N1(R,-1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2(R,-1,nu,vt,rho,z);
			terms += 3 * KerrAmendment::N3(R,-1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N4(R,-1,nu,vt,rho,z);
			break;
		}
		case 1: {
			terms += 3 * KerrAmendment::N1(R, 1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2(R, 1,nu,vt,rho,z);
			terms -= 3 * KerrAmendment::N3(R, 1,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N4(R, 1,nu,vt,rho,z);
			break; 
		}
		case -3: {
			terms += 1 * KerrAmendment::N1(R,-3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2(R,-3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N3(R,-3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N4(R,-3,nu,vt,rho,z);
			break; 
		}
		case 3: { 
			terms -= 1 * KerrAmendment::N1(R, 3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2(R, 3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N3(R, 3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N4(R, 3,nu,vt,rho,z);
			break; 
		}

		default: terms = 0;
	}
	
	return terms;
}

double KerrAmendment::N1 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);

	return 3 * m * jn(m, nu*rho) * i1 * i1 * i1_perp;
}

double KerrAmendment::N2 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	return m * jn(m, nu*rho) * (i2-i1) * (i1_perp*(i2-i1) + 2*i1*(i2_perp-i1_perp));
}

double KerrAmendment::N3 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);
	double bessel_diff = jn(m-1, rho*nu) - jn(m+1, rho*nu);

	return -1.5 * nu * rho * bessel_diff * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp);
}

double KerrAmendment::N4 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);
	double bessel_diff = jn(m-1,rho*nu) - jn(m+1,rho*nu);
	double imult_perp = i1 * (i2_perp - i1_perp) + 2 * i1_perp * (i2 - i1);

	return -0.5 * nu * rho * bessel_diff * i1 * imult_perp;
}

/* partial derivatives of I1 and I2 by time */

double KerrAmendment::int_bessel_011_perp (double vt, double z, double rho, double R)
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;
 
	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal!");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double first = (rho2 - R2) * (rho2 - R2) / vt_z;
	first /= std::sqrt( (rho+R)*(rho+R) - vt_z );
	first /= std::sqrt( vt_z - (rho-R)*(rho-R) );

	double second = 2 * (rho2 + R2) - vt_z;
	second /= std::sqrt(4*rho2*R2 - (vt_z - rho2 - R2)*(vt_z - rho2 - R2));

	return vt * (first - second) / (2*M_PI*rho2);
}

double KerrAmendment::int_bessel_001_perp (double vt, double z, double rho, double R)
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal!");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = - vt / vt_z;
	res *= vt_z - rho2 + R2;
	res /= std::sqrt(4*rho2*vt_z - (vt_z+rho2-R2) * (vt_z+rho2-R2) );
	return res / M_PI;
}
