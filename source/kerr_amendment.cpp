//
//  kerr_amendment.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

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

/* */

KerrAmendment::KerrAmendment (MissileField* field, KerrMedium* medium, UniformPlainDisk* source)
: NonlinearField(field, medium) 
{ 
	this->linear_field = field;
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
}

double KerrAmendment::electric_rho (double vt, double rho, double phi, double z) const
{
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-50) return 0;

	/* double eps0 = NonlinearMedium::EPS0;
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu0 = NonlinearMedium::MU0;
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = mu0 * mu_r / eps0 * eps_r;
	double A0 = this->A0;
	double coeff = kerr * A0 * A0 * A0 * em_relation * em_relation / 128; */

	/* Monte-Carlo Inegration 
	auto field = [this, rho, phi, z] (double vt) {

		auto mode = [this, vt, rho, phi, z] (int m, double varrho, double z_perp, double vt_perp, double nu) {
			double mu_r = nl_medium->relative_permeability(vt,z);
			double sqrt_eps0 = std::sqrt(NonlinearMedium::EPS0);
			double G = this->riemann(nu, vt_perp - vt, z_perp - z);
			double res = std::cos(m * phi);
			res *= m * G * jn(m, nu * rho);
			res /= rho * std::sqrt(nu);
			res *= varrho * this->im_modal_source(m,nu,vt,varrho,z);
			return res * mu_r / sqrt_eps0;
		};

		auto modes_sum = [this, vt, rho, phi, z, mode] (double varrho, double z_perp, double vt_perp, double nu) {
			double sum = 0;
			sum += mode(-1, varrho, z_perp, vt_perp, nu);
			sum += mode(1, varrho, z_perp, vt_perp, nu);
			sum += mode(-3, varrho, z_perp, vt_perp, nu);
			sum += mode(3, varrho, z_perp, vt_perp, nu);
			return sum;
		};

		std::vector<std::pair<double,double>> limits;
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		MonteCarlo integral = MonteCarlo(10e8, limits);
		return integral.value(modes_sum);
	};

	return coeff * Math::derivative(field,vt); */

	/* auto field = [this, rho, phi, z] (double vt) {

		auto mode = [this, vt, rho, phi, z] (int m, double rho_perp, double z_perp, double vt_perp, double nu) {
			double mu_r = nl_medium->relative_permeability(vt,z);
			double sqrt_eps0 = std::sqrt(NonlinearMedium::EPS0);
			double G = this->riemann(nu, vt_perp - vt, z_perp - z);
			double res = std::cos(m * phi);
			res *= m * G * jn(m, nu * rho);
			res /= rho * std::sqrt(nu);
			res *= rho_perp * this->im_modal_source(m,nu,vt,rho_perp,z);
			return res;
		};

		auto modes_sum = [mode] (double rho_perp, double z_perp, double vt_perp, double nu) {
			double sum = 0;
			sum += mode(-1, rho_perp, z_perp, vt_perp, nu);
			sum += mode( 1, rho_perp, z_perp, vt_perp, nu);
			sum += mode(-3, rho_perp, z_perp, vt_perp, nu);
			sum += mode( 3, rho_perp, z_perp, vt_perp, nu);
			return sum;
		};

		std::vector< std::tuple<double,std::size_t,double> > limits;
		limits.push_back(std::make_tuple(0, 4e3, 1e3)); // rho_perp
		limits.push_back(std::make_tuple(0, 40, 2));    // z_perp
		limits.push_back(std::make_tuple(0, 40, 2));    // vt_perp
		limits.push_back(std::make_tuple(0, 40, 2)); 	// nu_perp
		SimpsonMultiDim integral = SimpsonMultiDim(limits);
		return integral.value(modes_sum);
	};

	return coeff * Math::derivative(field,vt); */

	throw std::logic_error("KerrAmendment::electric_rho is not implemented");
}

double KerrAmendment::electric_phi (double vt, double rho, double phi, double z) const
{
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-50) return 0;
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::electric_phi is not implemented");
}

double KerrAmendment::electric_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

double KerrAmendment::magnetic_rho (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-50) return 0;
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-50) return 0;
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(phi); UNUSED(rho); UNUSED(z);
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-50) return 0;
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

double KerrAmendment::riemann (double nu, double vt_diff, double z_diff) const
{ 
	double distance = vt_diff * vt_diff - z_diff * z_diff;
	if (distance > 0) distance = std::sqrt(distance);
	else std::invalid_argument("Interval is not legal!");

	return j0(nu * distance);
}


double KerrAmendment::im_modal_source_sum (double nu, double vt, double rho, double z) const
{
	double terms = 0;

	terms -= 3 * KerrAmendment::N1(-1,nu,vt,rho,z);
	terms -= 1 * KerrAmendment::N2(-1,nu,vt,rho,z);
	terms += 3 * KerrAmendment::N4(-1,nu,vt,rho,z);
	terms += 1 * KerrAmendment::N5(-1,nu,vt,rho,z);

	terms -= 3 * KerrAmendment::N1( 1,nu,vt,rho,z);
	terms -= 1 * KerrAmendment::N2( 1,nu,vt,rho,z);
	terms -= 3 * KerrAmendment::N4( 1,nu,vt,rho,z);
	terms -= 1 * KerrAmendment::N5( 1,nu,vt,rho,z);

	terms -= 1 * KerrAmendment::N1(-3,nu,vt,rho,z);
	terms += 1 * KerrAmendment::N2(-3,nu,vt,rho,z);
	terms -= 1 * KerrAmendment::N4(-3,nu,vt,rho,z);
	terms += 1 * KerrAmendment::N5(-3,nu,vt,rho,z);

	terms += 1 * KerrAmendment::N1( 3,nu,vt,rho,z);
	terms += 1 * KerrAmendment::N2( 3,nu,vt,rho,z);
	terms -= 1 * KerrAmendment::N4( 3,nu,vt,rho,z);
	terms += 1 * KerrAmendment::N5( 3,nu,vt,rho,z);
	
	return terms;
}

double KerrAmendment::im_modal_source (int m, double nu, double vt, double rho, double z) const
{
	double terms = 0;

	switch (m) {
		case -1: { 
			terms -= 3 * KerrAmendment::N1(-1,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2(-1,nu,vt,rho,z);
			terms += 3 * KerrAmendment::N4(-1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5(-1,nu,vt,rho,z);
			break;
		}
		case 1: {
			terms -= 3 * KerrAmendment::N1( 1,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2( 1,nu,vt,rho,z);
			terms -= 3 * KerrAmendment::N4( 1,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N5( 1,nu,vt,rho,z);
			break; 
		}
		case -3: {
			terms -= 1 * KerrAmendment::N1(-3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2(-3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N4(-3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5(-3,nu,vt,rho,z);
			break; 
		}
		case 3: { 
			terms += 1 * KerrAmendment::N1( 3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2( 3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N4( 3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5( 3,nu,vt,rho,z);
			break; 
		}

		default: terms = 0;
	}
	
	return terms;
}

double KerrAmendment::N1 (int m, double nu, double vt, double rho, double z) const
{
	double vt_z = std::sqrt(vt * vt - z * z);
	double i1 = this->linear_field->int_bessel_011(vt_z, rho, this->R);
	double i1_perp = this->int_bessel_011_perp(vt, z, rho, this->R);

	return jn(m, nu*rho) * i1 * i1 * i1_perp;
}

double KerrAmendment::N2 (int m, double nu, double vt, double rho, double z) const
{
	double vt_z = std::sqrt(vt * vt - z * z);
	double i1 = this->linear_field->int_bessel_011(vt_z, rho, this->R);
	double i2 = this->linear_field->int_bessel_001(vt_z, rho, this->R);
	double i1_perp = this->int_bessel_011_perp(vt, z, rho, this->R);
	double i2_perp = this->int_bessel_001_perp(vt, z, rho, this->R);

	return (i2-i1) * (i1_perp*(i2-i1) + 2*i1*(i2_perp-i1_perp));
}

double KerrAmendment::N4 (int m, double nu, double vt, double rho, double z) const
{
	double vt_z = std::sqrt(vt * vt - z * z);
	double i1 = this->linear_field->int_bessel_011(vt_z, rho, this->R);
	double i2 = this->linear_field->int_bessel_001(vt_z, rho, this->R);
	double i1_perp = this->int_bessel_011_perp(vt, z, rho, this->R);
	double i2_perp = this->int_bessel_001_perp(vt, z, rho, this->R);
	double bessel_diff = jn(m-1, rho*nu) - jn(m+1, rho*nu);

	return (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp) * bessel_diff / 2;
}

double KerrAmendment::N5 (int m, double nu, double vt, double rho, double z) const
{
	double vt_z = std::sqrt(vt * vt - z * z);
	double i1 = this->linear_field->int_bessel_011(vt_z, rho, this->R);
	double i2 = this->linear_field->int_bessel_001(vt_z, rho, this->R);
	double i1_perp = this->int_bessel_011_perp(vt, z, rho, this->R);
	double i2_perp = this->int_bessel_001_perp(vt, z, rho, this->R);
	double bessel_diff = jn(m-1,rho*nu) - jn(m+1,rho*nu);
	double imult_perp = i1 * (i2_perp - i1_perp) + 2 * i1_perp * (i2 - i1);

	return imult_perp * bessel_diff / 2;
}

/* */

double KerrAmendment::int_bessel_011_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = this->R * this->R;

	if (vt_z <= 0) throw std::invalid_argument("ct-z <= 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (this->R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (this->R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (this->R > rho + std::sqrt(vt_z)) return 0.0;

	double res = vt / rho2;
	res *= (rho2 - R2) * (rho2 - R2) / vt_z;
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	return res / (2*M_PI);
}

double KerrAmendment::int_bessel_001_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = this->R * this->R;

	if (vt_z <= 0) throw std::invalid_argument("ct-z <= 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (this->R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (this->R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (this->R > rho + std::sqrt(vt_z)) return 0.0;

	double res = - vt / vt_z;
	res *= vt_z - rho2 + R2;
	res /= std::sqrt(4*rho2*vt_z - (vt_z+rho2-R2) * (vt_z+rho2-R2) );
	return res / M_PI;
}
