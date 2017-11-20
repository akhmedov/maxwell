//
//  uniform_disk_current.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "uniform_disk_current.hpp"

// Linear Homogeneous Medium

Homogeneous::Homogeneous(double realative_mu, double realative_eps)
{
	if (realative_mu > 0)
		this->mu_r = realative_mu;
	else
		std::invalid_argument("Relative permeability must be positive.");
	
	if (realative_eps > 0)
		this->eps_r = realative_eps;
	else
		std::invalid_argument("Relative permittivity must be positive.");
}

double Homogeneous::relative_permittivity (double ct, double z) const
{
	UNUSED(ct); UNUSED(z);
	return Homogeneous::eps_r;
}

double Homogeneous::relative_permeability (double ct, double z) const
{
	UNUSED(ct); UNUSED(z);
	return Homogeneous::mu_r;
}

// Linear umiform electric current

UniformPlainDisk::UniformPlainDisk(double disk_radius, double magnitude)
{
	if (disk_radius > 0) this->R = disk_radius;
	else std::invalid_argument("Relative permittivity must be positive.");
	
	if (magnitude > 0) this->A0 = magnitude;
	else std::invalid_argument("Electic current magnitude must be positive.");
}

double UniformPlainDisk::rho (double ct, double rho, double phi, double z) const
{
	double coeffisient = Math::kronecker_delta(z,0) * Math::heaviside_sfunc(ct,0);
	coeffisient *= Math::heaviside_sfunc(rho,0) - Math::heaviside_sfunc(rho,R);
	return coeffisient * std::cos(phi);
}

double UniformPlainDisk::phi (double ct, double rho, double phi, double z) const
{
	double coeffisient = Math::kronecker_delta(z,0) * Math::heaviside_sfunc(ct,0);
	coeffisient *= Math::heaviside_sfunc(rho,0) - Math::heaviside_sfunc(rho,R);
	return - coeffisient * std::sin(phi);
}

double UniformPlainDisk::z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

double UniformPlainDisk::get_magnitude() const
{
	return this->A0;
}

double UniformPlainDisk::get_disk_radius() const
{
	return this->R;
}

// Electric field

MissileField::MissileField (UniformPlainDisk* source, Homogeneous* medium)
: LinearField(source,medium)
{
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
}

/* MissileField::MissileField (UniformPlainDisk* source, KerrMedium* medium)
: LinearField(source,medium)
{
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
} */

void MissileField::set_yterms_num (std::size_t number)
{
	this->STATIC_TERMS_NUMBER = number;
}

std::size_t MissileField::get_yterms_num () const
{
	return this->STATIC_TERMS_NUMBER;
}

double MissileField::electric_rho (double ct, double rho, double phi, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	
	double value = this->A0 / 2;
	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(ct,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(ct,z));
	value *= this->int_bessel_011( std::sqrt(vt_z.get_d()), rho, this->R) * std::cos(phi);
	return value;
}

double MissileField::electric_phi (double ct, double rho, double phi, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), rho, this->R);
	double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), rho, this->R);
	double value = this->A0 / 2;
	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(ct,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(ct,z));
	value *= - (i2 - i1) * std::sin(phi);
	return value;
}

double MissileField::electric_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

// Magnetic field

double MissileField::magnetic_rho (double ct, double rho, double phi, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	mpf_class ct_z = mpf_class(ct * ct);
	mpf_div(ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (ct_z.get_d() <= 0) return 0;

	double value = this->A0 / 2;
	value *= (this->int_lommel_001(ct, rho, z) - this->int_lommel_011(ct, rho, z)) * std::sin(phi);
	return value;
}

double MissileField::magnetic_phi (double ct, double rho, double phi, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	mpf_class ct_z = mpf_class(ct * ct);
	mpf_div(ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (ct_z.get_d() <= 0) return 0;

	double value = this->A0 / 2;
	value *= this->int_lommel_011(ct, rho, z) * std::cos(phi);
	return value;
}

double MissileField::magnetic_z (double ct, double rho, double phi, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	mpf_class ct_z = mpf_class(ct * ct);
	mpf_div(ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (ct_z.get_d() <= 0) return 0;

	double value = this->A0;
	value *= this->int_lommel_111(ct, rho, z) * std::sin(phi);
	return value;
}

double MissileField::static_magnitude (double z) const
{
	double R2 = this->R * this->R;
	double z2 = z*z;
	double res = - this->A0;
	double numer = std::sqrt(R2 - z2) - z;
	double denumer = std::sqrt(R2 - z2) + z;
	double sum = 0;
	for (std::size_t m = 0; m < STATIC_TERMS_NUMBER; m++)
		sum += std::pow(-numer/denumer,m+1);
	return res * sum;
}

double MissileField::static_magnitude (double rho, double phi, double z, double eps) const
{
	if ((z + rho) < 0.05) return this->A0/4;
	double R2 = this->R * this->R;
	double z2 = z*z;
	double ct = rho + std::sqrt(R2 + z2) + eps;
	return this->magnetic_y(ct,rho,phi,z);
}

// integrals

double MissileField::int_bessel_001 (double sqrt_vt_z, double rho, double R)
{
	if (sqrt_vt_z == 0) throw std::invalid_argument("ct-z = 0 is not allowed");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (rho == 0) {
		if (R < sqrt_vt_z) return 0;
		if (R == sqrt_vt_z) return 0.5;
		if (R > sqrt_vt_z) return 1;
	}

	if (R < std::abs(rho - sqrt_vt_z)) return 0.0;
	if (R > rho + sqrt_vt_z) return 1.0;

	double R2 = R * R;
	double rho2 = rho * rho;
	double vt_z = sqrt_vt_z * sqrt_vt_z;
	double numer = vt_z + rho2 - R2;
	double denumer = 2 * rho * sqrt_vt_z;
	return std::acos(numer/denumer) / M_PI;
}

double MissileField::int_bessel_011 (double sqrt_vt_z, double rho, double R)
{	
	if (sqrt_vt_z == 0) throw std::invalid_argument("ct-z = 0 is not allowed");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");
	
	if (rho == 0) {
		if (R < sqrt_vt_z) return 0;
		if (R == sqrt_vt_z) return 0.25;
		if (R > sqrt_vt_z) return 0.5;
	}
	
	if (R < std::abs(rho - sqrt_vt_z)) return 0;
	if (R > rho + sqrt_vt_z) return 0.5;
	
	double res = 0;
	double R2 = R * R;
	double rho2 = rho * rho;
	double vt_z = sqrt_vt_z * sqrt_vt_z;
	double frac = (rho - R) * (rho - R) / (rho + R) / (rho + R);
	double numer = (rho + R) * (rho + R) - vt_z;
	double denumer = vt_z  - (rho - R) * (rho - R);
	res += (rho2 + R2) * std::acos((vt_z - rho2 - R2)/ (2 * rho * R));
	res -= std::sqrt(4*rho2*R2 - (rho2 + R2 - vt_z) * (rho2 + R2 - vt_z));
	res -= 2 * std::abs(rho2 - R2) * std::atan(std::sqrt(frac*numer/denumer));
	return res / (4 * M_PI * rho2);
}

double MissileField::int_lommel_001 (double ct, double rho, double z) const
{
	double eps_r = medium->relative_permittivity(ct,z);
	double mu_r = medium->relative_permeability(ct,z);
	double R2 = this->R * this->R;

	mp_bitcnt_t bitrate = mp_bitcnt_t(256);
	mpf_set_default_prec(bitrate);

	mpf_class ctpz, ctmz, ctz;
	mpf_div(ctpz.get_mpf_t(), mpf_class(ct).get_mpf_t(), mpf_class(std::sqrt(eps_r * mu_r)).get_mpf_t());
	mpf_add(ctpz.get_mpf_t(), ctpz.get_mpf_t(), mpf_class(z).get_mpf_t());
	mpf_div(ctmz.get_mpf_t(), mpf_class(ct).get_mpf_t(), mpf_class(std::sqrt(eps_r * mu_r)).get_mpf_t());
	mpf_add(ctmz.get_mpf_t(), ctmz.get_mpf_t(), mpf_class(-z).get_mpf_t());
	mpf_div(ctz.get_mpf_t(), ctmz.get_mpf_t(), ctpz.get_mpf_t());
	mpf_class ctz_in_pow(ctz);

	mpf_class ct_z = mpf_class(ct * ct);
	mpf_div(ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	
	mpf_class yacob_arg, yacob;

	if (rho == 0) {
		if (R < std::sqrt(ct_z.get_d())) {
			mpf_class sum = mpf_class(0);
			for (std::size_t m = 0; m <= this->STATIC_TERMS_NUMBER; m++) {
				yacob = Math::yacobi_polinom(m, 1 - 2 * R2 / ct_z);
				mpf_mul(yacob.get_mpf_t(), yacob.get_mpf_t(), ctz_in_pow.get_mpf_t());
				mpf_add(sum.get_mpf_t(), sum.get_mpf_t(), yacob.get_mpf_t());
				mpf_mul(ctz_in_pow.get_mpf_t(), ctz_in_pow.get_mpf_t(), ctz.get_mpf_t());
			}
			mpf_class res = mpf_class(2 * R2);
			mpf_mul(res.get_mpf_t(), res.get_mpf_t(), sum.get_mpf_t());
			mpf_div(res.get_mpf_t(), res.get_mpf_t(), ct_z.get_mpf_t());
			return res.get_d();
		}
		if (R >= std::sqrt(ct_z.get_d())) return 1;
	}

	throw std::logic_error("Not implemented for non zero rho!");
}

double MissileField::int_lommel_011 (double ct, double rho, double z) const
{
	return this->int_lommel_001(ct, rho, z) / 2;
}

double MissileField::int_lommel_111 (double ct, double rho, double z) const
{
	UNUSED(ct); UNUSED(rho); UNUSED(z);
	throw std::logic_error("Not implemented");
}
