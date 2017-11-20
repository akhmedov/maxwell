//
//  kerr_amendment.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

KerrMedium::KerrMedium (double realative_mu, double realative_eps, double chi3_electric, double conduct) 
: Homogeneous(realative_mu, realative_eps) {
	if (chi3_electric > 0)
		this->chi3 = chi3_electric;
	else
		std::invalid_argument("Kerr coefficient must be positive.");

	
	if (conduct > 0) 
		this->sigma = conduct;
	else
		std::invalid_argument("Kerr coefficient must be positive.");
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
		case 3: return this->chi3;
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

KerrAmendment::KerrAmendment (MissileField* field, KerrMedium* medium)
: MissileField(*field), NonlinearField(field, medium) { }

double KerrAmendment::electric_rho (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::electric_phi is not implemented");
}

double KerrAmendment::electric_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::electric_phi is not implemented");
}

double KerrAmendment::electric_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

double KerrAmendment::magnetic_rho (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

std::complex<double> KerrAmendment::complex_electric_rho (double vt, double rho, double phi, double z) const
{

	auto mode = [this, vt, rho, phi, z] (std::size_t m, double nu) {
		std::complex<double> res =  - std::cos(m * phi) - 1i * std::sin(m * phi);
		res /= std::sqrt( LinearMedium::EPS0 );
		res *= jn(m, nu * rho) * this->evolution_Vh(m, nu, vt, z).imag();
		res /= rho * std::sqrt(nu);
		return res;
	};

	auto re_modes = [mode] (double nu) { 
		return mode(-1, nu).real() + mode(1, nu).real() + 
		mode(3, nu).real() + mode(-3, nu).real(); 
	};

	auto im_modes = [mode] (double nu) { 
		return mode(-1, nu).imag() + mode(1, nu).imag() + 
		mode(3, nu).imag() + mode(-3, nu).imag(); 
	};

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	return I.value(0, bais, re_modes) + I.value(0, bais, im_modes);
}

std::complex<double> KerrAmendment::complex_electric_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

std::complex<double> KerrAmendment::complex_electric_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::electric_z is not implemented");
}

std::complex<double> KerrAmendment::complex_magnetic_rho (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

std::complex<double> KerrAmendment::complex_magnetic_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

std::complex<double> KerrAmendment::complex_magnetic_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

std::complex<double> KerrAmendment::evolution_Ih (long m, double nu, double vt, double z) const
{
	auto Re_hm = [&] (double val_z) {
		return this->evolution_h(m, nu, vt, val_z).real();
	};

	auto Im_hm = [&] (double val_z) {
		return this->evolution_h(m, nu, vt, val_z).imag();
	};

	return Math::derivative(Re_hm, z) + 1i * Math::derivative(Im_hm, z);
}

std::complex<double> KerrAmendment::evolution_Vh (long m, double nu, double vt, double z) const
{
	double mu_r = nl_medium->relative_permeability(vt,z);

	auto Im_hm = [&] (double val_vt) {
		return this->evolution_h(m, nu, val_vt, z).imag();
	};

	auto Re_hm = [&] (double val_vt) {
		return this->evolution_h(m, nu, val_vt, z).real();
	};

	return - mu_r * ( Math::derivative(Re_hm, vt) + 1i * Math::derivative(Im_hm, vt) );
}

std::complex<double> KerrAmendment::evolution_h (long m, double nu, double vt, double z) const
{
	size_t bais = 10e5;
	Simpson It = Simpson(10*bais);

	auto func1 = [&] (double vt_perp) {
		std::cout << nu << ' ' << vt_perp << std::endl;
		Simpson Iz = Simpson(10*bais);
		auto func2 = [&] (double z_perp) {
			if ( vt_perp - z_perp < 10e-10) return 0.0;
			double G = this->riemann(nu, vt - vt_perp, z - z_perp);
			// if (G < 10e-10) return 0.0;
			std::complex<double> jm = this->modal_source(m, nu, vt_perp, z_perp);
			return jm.imag() * G;
		};
		return Iz.value(vt_perp, vt_perp + 2 * this->R, func2);
	};

	return It.value(10e-10, bais, func1);
}

double KerrAmendment::riemann (double nu, double vt_diff, double z_diff) const
{ 
	double distance = vt_diff * vt_diff - z_diff * z_diff;
	if (distance > 0) distance = std::sqrt(distance);
	else std::invalid_argument("Interval is not legal!");

	double eps_r = nl_medium->relative_permittivity(vt_diff,z_diff);
	double mu_r = nl_medium->relative_permeability(vt_diff,z_diff);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);

	double interval = vt_diff - z_diff;
	if (interval > 0) return j0(nu * distance) * velocity / 2;
	// else if (interval == 0) return velocity / 4; 
	else return 0;
}

std::complex<double> KerrAmendment::modal_source (long m, double nu, double ct, double z) const
{
	double terms = 0;

	switch (m) {
		case -1: { 
			terms -= 3 * KerrAmendment::N1(-1,nu,ct,z);
			terms -= 1 * KerrAmendment::N2(-1,nu,ct,z);
			terms -= 1 * KerrAmendment::N3(-1,nu,ct,z);
			terms += 3 * KerrAmendment::N4(-1,nu,ct,z);
			terms += 1 * KerrAmendment::N5(-1,nu,ct,z);
			terms += 1 * KerrAmendment::N6(-1,nu,ct,z);
			break;
		}
		case 1: {
			terms -= 3 * KerrAmendment::N1(1,nu,ct,z);
			terms -= 1 * KerrAmendment::N2(1,nu,ct,z);
			terms -= 1 * KerrAmendment::N3(1,nu,ct,z);
			terms -= 3 * KerrAmendment::N4(1,nu,ct,z);
			terms -= 1 * KerrAmendment::N5(1,nu,ct,z);
			terms -= 1 * KerrAmendment::N6(1,nu,ct,z);
			break; 
		}
		case -3: {
			terms -= 1 * KerrAmendment::N1(-3,nu,ct,z);
			terms += 1 * KerrAmendment::N2(-3,nu,ct,z);
			terms -= 1 * KerrAmendment::N4(-3,nu,ct,z);
			terms += 1 * KerrAmendment::N5(-3,nu,ct,z);
			break; 
		}
		case 3: { 
			terms += 1 * KerrAmendment::N1(3,nu,ct,z);
			terms += 1 * KerrAmendment::N2(3,nu,ct,z);
			terms -= 1 * KerrAmendment::N4(3,nu,ct,z);
			terms += 1 * KerrAmendment::N5(3,nu,ct,z);
			break; 
		}

		default: terms = 0;
	}
	
	return 1i * terms;
}

double KerrAmendment::N1 (long m, double nu, double ct, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	double vt = ct / std::sqrt(eps_r * mu_r);

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, rho, this->R);
		return jn(m,rho*nu) * i1 * i1 * i1_perp;
	};

	double res = 3 * m * A3 * chi * std::sqrt(NonlinearMedium::MU0);
	res /= 64 * std::sqrt(nu);
	res *= em_relation * em_relation * em_relation;
	res *= I.value(0, bais, func);
	return res;
}

double KerrAmendment::N2 (long m, double nu, double ct, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), rho, this->R);
		double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, rho, this->R);
		double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, rho, this->R);
		double res = i1_perp * (i2 - i1) + 2 * i1 * (i2_perp - i1_perp);
		return (i2 - i1) * res;
	};

	double res = m * A3 * chi * std::sqrt(NonlinearMedium::MU0);
	res /= 64 * std::sqrt(nu);
	res *= em_relation * em_relation * em_relation;
	res *= I.value(0, bais, func);
	return res;
}

double KerrAmendment::N3 (long m, double nu, double ct, double z) const
{
	double sigma = nl_medium->conductivity(ct, z);
	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		return i1 * jn(m,rho*nu);
	};

	double res = m * this->A0 * sigma * std::sqrt(NonlinearMedium::MU0);
	res /= 4 * std::sqrt(nu);
	res *= em_relation * I.value(0, bais, func);
	return res;
}

double KerrAmendment::N4 (long m, double nu, double ct, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), rho, this->R);
		double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, rho, this->R);
		double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, rho, this->R);
		double bessel = jn(m-1,rho*nu) - jn(m+1,rho*nu);
		return rho * (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp) * bessel;
	};

	double res = A3 * chi * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 128;
	res *= em_relation * em_relation * em_relation;
	res *= I.value(0, bais, func);
	return res;
}

double KerrAmendment::N5 (long m, double nu, double ct, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), rho, this->R);
		double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, rho, this->R);
		double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, rho, this->R);
		double bessel = jn(m-1,rho*nu) - jn(m+1,rho*nu);
		double mult = i1 * (i2_perp - i1_perp) + 2 * i1_perp * (i2 - i1);
		return rho * mult * bessel;
	};

	double res = 3 * A3 * chi * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 128;
	res *= em_relation * em_relation * em_relation;
	res *= I.value(0, bais, func);
	return res;
}

double KerrAmendment::N6 (long m, double nu, double ct, double z) const
{
	double sigma = nl_medium->conductivity(ct, z);
	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	auto func = [&] (double rho) { 
		double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
		double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), rho, this->R);
		double bessel = jn(m-1,rho*nu) - jn(m+1,rho*nu);
		return rho * (i2 - i1) * bessel;
	};

	double res = this->A0 * sigma * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 8;
	res *= em_relation * I.value(0, bais, func);
	return res;
}

/* */

double KerrAmendment::vint_bessel_011_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal");
	if (vt_z == 0) throw std::invalid_argument("ct-z = 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = vt / rho2;
	res *= (rho2 - R2) * (rho2 - R2) / vt_z;
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	return res / (2*M_PI);
}

double KerrAmendment::vint_bessel_001_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal");
	if (vt_z == 0) throw std::invalid_argument("ct-z = 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = - vt / vt_z;
	res *= vt_z - rho2 + R2;
	res /= std::sqrt(4*rho2*vt_z - (vt_z+rho2-R2) * (vt_z+rho2-R2) );
	return res / M_PI;
}
