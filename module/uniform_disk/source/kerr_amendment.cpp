//
//  kerr_amendment.cpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

KerrAmendment::KerrAmendment (double R, double A0, double eps_r, double mu_r, double chi3, Logger* log)
: TransientResponse(R,A0,eps_r,mu_r,log), kerr(chi3) {}

// ===============================================================================================

double KerrAmendment::current_x (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	const double current_rho = this->current_rho(event);
	const double current_phi = this->current_phi(event);
	return current_rho * std::cos(event.phi()) - current_phi * std::sin(event.phi());
}

double KerrAmendment::current_rho (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double em_relation = std::sqrt(MU0 * MU / EPS0 * EPS);
	double A3 = this->A0 * this->A0 * this->A0;

	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	double mult1 = i1_perp * (i2 - i1);
	double mult2 = i1 * (i2_perp - i1_perp);

	double coeff = A3 * em_relation * em_relation * em_relation / 8;
	double cos3sin0 = 3 * i1 * i1 * i1_perp * std::cos(phi) * std::cos(phi) * std::cos(phi);
	double cos1sin2 = (i2 - i1) * (mult1 + 2 * mult2) * std::cos(phi) * std::sin(phi) * std::sin(phi);
	return coeff * (cos3sin0 + cos1sin2);
}

double KerrAmendment::current_phi (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double em_relation = std::sqrt(MU0 * MU / EPS0 * EPS);
	double A3 = this->A0 * this->A0 * this->A0;

	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	double mult1 = i1_perp * (i2 - i1);
	double mult2 = i1 * (i2_perp - i1_perp);

	double coeff = A3 * em_relation * em_relation * em_relation / 8;
	double cos0sin3 = 3 * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp) * std::sin(phi) * std::sin(phi) * std::sin(phi);
	double cos2sin1 = i1 * (2 * mult1 + mult2) * std::cos(phi) * std::cos(phi) * std::sin(phi);
	return - coeff * (cos0sin3 + cos2sin1);
}

// ===============================================================================================

double KerrAmendment::observed_from (const Point::Cylindrical&) const
{
	throw std::logic_error("KerrAmendment::observed_from is not implemented!");
}

double KerrAmendment::observed_to (const Point::Cylindrical&) const
{
	throw std::logic_error("KerrAmendment::observed_to is not implemented!");
}

double KerrAmendment::electric_x (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double vt_z = vt - z;
	if (z == 0) return 0;
	if (vt_z < 1e-9) return 0;

	double A3 = this->A0 * this->A0 * this->A0;
	double em_relation = MU0 * MU / EPS0 * EPS;
	double coeff = this->kerr * A3 * em_relation * em_relation / 128;

	double integral;

	try {

		integral = SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, 2*R,
			[this, vt, rho, phi, z] (double z_perp) {

				try {

					return SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, 2*R,
						[this, vt, rho, phi, z, z_perp] (double rho_perp) {

							double max_vt = vt - z + z_perp; // grater then zero
							double max_nu = PERIODS_NU * std::abs(rho - rho_perp);
							double period_nu =  NODES_NU * (rho + rho_perp);
							std::size_t nodes_nu = max_nu * period_nu;
							if (!nodes_nu) nodes_nu = 1;

							double delta_sum;

							try  {

								delta_sum = SimpsonRunge(nodes_nu, MAX_ERROR).value(0, max_nu,
									[this, rho, phi, max_vt, rho_perp, z_perp] (double nu) {

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

								if (this->global_log) {
									std::string mesg = KerrAmendment::INTEGRAL_WARNING;
									Point::SpaceTime<Point::Cylindrical> event{vt,rho,phi,z};
									mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "Z'");
									mesg = std::regex_replace(mesg, std::regex("\\$POINT"), event.to_str());
									this->global_log->warning(mesg);
								}

								delta_sum = not_trusted;
							}
							
							if (std::isnan(delta_sum)) delta_sum = 0;

							double step_sum = 0;
							try {

								step_sum = SimpsonRunge(MIN_NODES, MAX_ERROR).value(0, max_vt,
									[this, vt, rho, phi, z, rho_perp, z_perp] (double vt_perp) {

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
											[this, rho, phi, rho_perp, z_perp, vt_perp, casual, delta_vt] (double nu) {
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

										if (this->global_log) {
											std::string mesg = KerrAmendment::INTEGRAL_WARNING;
											Point::SpaceTime<Point::Cylindrical> event{vt,rho,phi,z};
											mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "Z'");
											mesg = std::regex_replace(mesg, std::regex("\\$POINT"), event.to_str());
											this->global_log->warning(mesg);
										}

										res = not_trusted;
									}

									if (std::isnan(res)) res = 0;
									return res;

								} ); // <---------------------------------------------------------------------------------------- d tau'

							} catch (double not_trusted) {

								if (this->global_log) {
									std::string mesg = KerrAmendment::INTEGRAL_WARNING;
									Point::SpaceTime<Point::Cylindrical> event{vt,rho,phi,z};
									mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "Z'");
									mesg = std::regex_replace(mesg, std::regex("\\$POINT"), event.to_str());
									this->global_log->warning(mesg);
								}

								step_sum = not_trusted;
							}

						return delta_sum - step_sum;

					} ); //	<------------------------------------------------------------------------------------------------------- d rho'

				} catch (double not_trusted) {

					if (this->global_log) {
						std::string mesg = KerrAmendment::INTEGRAL_WARNING;
						Point::SpaceTime<Point::Cylindrical> event{vt,rho,phi,z};
						mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "Z'");
						mesg = std::regex_replace(mesg, std::regex("\\$POINT"), event.to_str());
						this->global_log->warning(mesg);
					}

					return not_trusted;
				}
		} );	// <---------------------------------------------------------------------------------------------------------------- d z'

	} catch (double not_trusted) {

		if (this->global_log) {
			std::string mesg = KerrAmendment::INTEGRAL_WARNING;
			Point::SpaceTime<Point::Cylindrical> tmp{vt,rho,phi,z};
			mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "Z'");
			mesg = std::regex_replace(mesg, std::regex("\\$POINT"), tmp.to_str());
			this->global_log->warning(mesg);
		}

		integral = not_trusted;
	}

	return - coeff * integral;
}

double KerrAmendment::electric_rho (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::electric_phi (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::electric_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	return 0;
}

double KerrAmendment::magnetic_rho (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
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
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);

	return 3 * m * jn(m, nu*rho) * i1 * i1 * i1_perp;
}

double KerrAmendment::N2 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);

	return m * jn(m, nu*rho) * (i2-i1) * (i1_perp*(i2-i1) + 2*i1*(i2_perp-i1_perp));
}

double KerrAmendment::N3 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);
	double bessel_diff = jn(m-1, rho*nu) - jn(m+1, rho*nu);

	return -1.5 * nu * rho * bessel_diff * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp);
}

double KerrAmendment::N4 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
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
