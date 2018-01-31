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
	if (kerr < 1e-10) return 0;

	double eps0 = NonlinearMedium::EPS0;
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu0 = NonlinearMedium::MU0;
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = mu0 * mu_r / eps0 * eps_r;
	double A0 = this->A0;
	double R = this->R;
	double coeff = kerr * A0 * A0 * A0 * em_relation * em_relation / 128;

	/* next section is not correct for arbirtaty inhomogenious medium */ 
	mpf_class vt_z = mpf_class(vt * vt);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	auto field = [R, rho, phi, z] (double vt) {

		auto mode = [R, vt, rho, phi, z] (int m, double rho_perp, double z_perp, double vt_perp, double nu) {
			double res = m * std::cos(m * phi);
			res *= rho ? jn(m, nu * rho) / rho * nu : 0.5;
			res *= KerrAmendment::riemann(nu, vt_perp - vt, z_perp - z);
			res *= KerrAmendment::im_modal_source(R,m,nu,vt,rho_perp,z);
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

	return coeff * Math::derivat3(field,vt);
}

double KerrAmendment::electric_phi (double vt, double rho, double phi, double z) const
{
	double kerr = this->nl_medium->relative_permittivity(vt,z,3);
	if (kerr < 1e-10) return 0;

	double eps0 = NonlinearMedium::EPS0;
	double eps_r = this->nl_medium->relative_permittivity(vt,z);
	double mu0 = NonlinearMedium::MU0;
	double mu_r = this->nl_medium->relative_permeability(vt,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = mu0 * mu_r / eps0 * eps_r;
	double A0 = this->A0;
	double R = this->R;
	double coeff = kerr * A0 * A0 * A0 * em_relation * em_relation / 128;

	/* next section is not correct for arbirtaty inhomogenious medium */ 
	mpf_class vt_z = mpf_class(vt * vt);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	auto field = [R, rho, phi, z] (double vt) {

		auto mode = [R, vt, rho, phi, z] (int m, double rho_perp, double z_perp, double vt_perp, double nu) {
			/* double res = m * std::cos(m * phi);
			res *= jn(m, nu * rho) * KerrAmendment::riemann(nu, vt_perp - vt, z_perp - z);
			res /= rho * std::sqrt(nu);
			res *= rho_perp * KerrAmendment::im_modal_source(R,m,nu,vt,rho_perp,z);
			return res; */
			throw std::logic_error("KerrAmendment::electric_phi is not implemented");
			return 0.0;
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

	return coeff * Math::derivat3(field,vt);
	
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

double KerrAmendment::riemann (double nu, double vt_diff, double z_diff)
{ 
	double distance = vt_diff * vt_diff - z_diff * z_diff;
	if (distance > 0) distance = std::sqrt(distance);
	else std::invalid_argument("Interval is not legal!");

	return j0(nu * distance);
}




double KerrAmendment::im_modal_source_sum (double R, double nu, double vt, double rho, double z)
{
	double mode_m3 = KerrAmendment::im_modal_source(R,-3,nu,vt,rho,z);
	double mode_m1 = KerrAmendment::im_modal_source(R,-1,nu,vt,rho,z);
	double mode_p1 = KerrAmendment::im_modal_source(R, 1,nu,vt,rho,z);
	double mode_p3 = KerrAmendment::im_modal_source(R, 3,nu,vt,rho,z);
	
	return mode_m3 + mode_m1 + mode_p1 + mode_p3;
}

double KerrAmendment::im_modal_source (double R, int m, double nu, double vt, double rho, double z)
{
	double terms = 0;

	switch (m) {
		case -1: { 
			terms += 3 * KerrAmendment::N1(R,-1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2(R,-1,nu,vt,rho,z);
			terms += 3 * KerrAmendment::N4(R,-1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5(R,-1,nu,vt,rho,z);
			break;
		}
		case 1: {
			terms += 3 * KerrAmendment::N1(R, 1,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N2(R, 1,nu,vt,rho,z);
			terms -= 3 * KerrAmendment::N4(R, 1,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N5(R, 1,nu,vt,rho,z);
			break; 
		}
		case -3: {
			terms += 1 * KerrAmendment::N1(R,-3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2(R,-3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N4(R,-3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5(R,-3,nu,vt,rho,z);
			break; 
		}
		case 3: { 
			terms -= 1 * KerrAmendment::N1(R, 3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N2(R, 3,nu,vt,rho,z);
			terms -= 1 * KerrAmendment::N4(R, 3,nu,vt,rho,z);
			terms += 1 * KerrAmendment::N5(R, 3,nu,vt,rho,z);
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

	return -m * (i2-i1) * (i1_perp*(i2-i1) + 2*i1*(i2_perp-i1_perp));
}

double KerrAmendment::N4 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);
	double bessel_diff = jn(m-1, rho*nu) - jn(m+1, rho*nu);

	return -1.5 * nu * rho * bessel_diff * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp);
}

double KerrAmendment::N5 (double R, int m, double nu, double vt, double rho, double z)
{
	double sqrt_vt_z = std::sqrt(vt * vt - z * z);
	double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(vt, z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(vt, z, rho, R);
	double bessel_diff = jn(m-1,rho*nu) - jn(m+1,rho*nu);
	double imult_perp = i1 * (i2_perp - i1_perp) + 2 * i1_perp * (i2 - i1);

	return -1.5 * nu * rho * bessel_diff * imult_perp;
}

/* partial derivatives of I1 and I2 by time */

double KerrAmendment::int_bessel_011_perp (double vt, double z, double rho, double R)
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z <= 0) throw std::invalid_argument("ct-z <= 0 is not legal");
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

	if (vt_z <= 0) throw std::invalid_argument("ct-z <= 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = - vt / vt_z;
	res *= vt_z - rho2 + R2;
	res /= std::sqrt(4*rho2*vt_z - (vt_z+rho2-R2) * (vt_z+rho2-R2) );
	return res / M_PI;
}
