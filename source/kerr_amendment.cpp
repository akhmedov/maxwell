//
//  kerr_amendment.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

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

KerrAmendment::KerrAmendment (MissileField* field, KerrMedium* medium, UniformPlainDisk* source)
: NonlinearField(field, medium) 
{ 
	this->linear_field = field;
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
}

double KerrAmendment::electric_x (double vt, double rho, double phi, double z) const
{
	double vt_z = vt - z;
	if (z == 0) return 0;
	if (vt_z < 1e-9) return 0;

	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double em_relation = NonlinearMedium::MU0 * mu_r / NonlinearMedium::EPS0 * eps_r;
	double coeff = kerr * this->A0 * this->A0 * this->A0 * em_relation * em_relation / 128;
	double R = this->R; 

	auto mode_sum = [R, vt, rho, phi, z] (double rho_perp, double z_perp) {
		
		double max_vt = vt - z + z_perp; // grater then zero

		auto delta_sum = [R, vt, rho, phi, z, max_vt, rho_perp, z_perp] (double nu) {
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
		};

		auto step_sum = [R, vt, rho, phi, z, rho_perp, z_perp] (double nu, double vt_perp) {
			double delta_vt = vt - vt_perp;
			double delta_z = z - z_perp;
			if (delta_vt - delta_z <= 0) return 0.0;
			double casual = std::sqrt(delta_vt*delta_vt - delta_z*delta_z);
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
		};

		double max_nu1 = 10 * std::abs(rho - rho_perp);
		double step_pints_nu = 7 * (rho + rho_perp) + 2;
		Simpson nu_int = Simpson(step_pints_nu * max_nu1);

		Simpson2D_line tau_nu_int = Simpson2D_line();
		tau_nu_int.first_limit(0, 50, max_vt);
		auto min_nu2 = [] (double vt_perp) { 
			return 0; 
		};
		auto max_nu2 = [rho, rho_perp, vt, z, z_perp] (double vt_perp) { 
			return 10 * std::abs(rho - rho_perp - std::sqrt(vt - vt_perp - z + z_perp)); 
		};
		tau_nu_int.second_limit(min_nu2, 700, max_nu2);
		
		double first = nu_int.value(0, max_nu1, delta_sum);
		if (std::isnan(first)) first = 0;
		double second = tau_nu_int.value(step_sum);
		if (std::isnan(second)) second = 0;
		return first - second;
	};

	std::vector< std::tuple<double,std::size_t,double> > limits;
	limits.push_back(std::make_tuple(0, 50, 2*R));
	limits.push_back(std::make_tuple(0, 50, 2*R));
	Simpson2D multi = Simpson2D(limits);
	return - coeff * multi.value(mode_sum);
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
