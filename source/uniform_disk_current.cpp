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
: LinearCurrent() {
	if (disk_radius > 0) this->R = disk_radius;
	else std::invalid_argument("Relative permittivity must be positive.");
	
	if (magnitude > 0) this->A0 = magnitude;
	else std::invalid_argument("Electic current magnitude must be positive.");
}

double UniformPlainDisk::time_shape (double vt) const
{
	return Math::heaviside_sfunc(vt,0);
}

double UniformPlainDisk::rho (double rho, double phi, double z) const
{
	double coeffisient = Math::kronecker_delta(z,0);
	coeffisient *= Math::heaviside_sfunc(rho,0) - Math::heaviside_sfunc(rho,R);
	return coeffisient * std::cos(phi);
}

double UniformPlainDisk::phi (double rho, double phi, double z) const
{
	double coeffisient = Math::kronecker_delta(z,0);
	coeffisient *= Math::heaviside_sfunc(rho,0) - Math::heaviside_sfunc(rho,R);
	return - coeffisient * std::sin(phi);
}

double UniformPlainDisk::z (double rho, double phi, double z) const
{
	UNUSED(phi); UNUSED(rho); UNUSED(z);
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

std::size_t MissileField::STATIC_TERMS_NUMBER = 100;

MissileField::MissileField (UniformPlainDisk* source, Homogeneous* medium)
: LinearField(source,medium)
{
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
}

void MissileField::set_yterms_num (std::size_t number)
{
	MissileField::STATIC_TERMS_NUMBER = number;
}

std::size_t MissileField::get_yterms_num ()
{
	return MissileField::STATIC_TERMS_NUMBER;
}

double MissileField::electric_rho (double ct, double rho, double phi, double z) const
{
	double vt_z = ct*ct - z*z;
	if (vt_z < 1e-9) return 0;
	double sqrt_vt_z = std::sqrt(vt_z);
	double R = this->R;
	double value = this->A0 / 2;
	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(ct,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(ct,z));

	#ifdef NUMERIC_PDISK_LINEAR_INT

		std::size_t bais = 1e5;
		Simpson I = Simpson(10*bais);
		
		auto ui_011 = [sqrt_vt_z, rho, R] (double nu) {
			if (nu == 0) return 0.0;
			if (rho == 0) return R * j0(nu*sqrt_vt_z) * j1(nu*R) / 2;
			return R * j0(nu*sqrt_vt_z) * j1(nu*rho) * j1(nu*R) / rho / nu;
		};

		value *= I.value(0, bais, ui_011) * std::cos(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */

		double i1 = this->int_bessel_011(sqrt_vt_z, rho, R);
		value *= i1 * std::cos(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double MissileField::electric_phi (double ct, double rho, double phi, double z) const
{
	double vt_z = ct*ct - z*z;
	if (vt_z < 1e-9) return 0;
	double sqrt_vt_z = std::sqrt(vt_z);
	double R = this->R;
	double value = this->A0 / 2;
	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(ct,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(ct,z));

	#ifdef NUMERIC_PDISK_LINEAR_INT

		std::size_t bais = 1e5;
		Simpson I = Simpson(10*bais);

		auto ui_001 = [sqrt_vt_z, rho, R] (double nu) { 
			return j0(nu*sqrt_vt_z) * j0(nu*rho) * j1(nu*R); 
		};
		
		auto ui_011 = [sqrt_vt_z, rho, R] (double nu) {
			if (nu == 0) return 0.0;
			if (rho == 0) return R * j0(nu*sqrt_vt_z) * j1(nu*R) / 2;
			return R * j0(nu*sqrt_vt_z) * j1(nu*rho) * j1(nu*R) / rho / nu;
		};

		value *= - (I.value(0, bais, ui_001) - I.value(0, bais, ui_011)) * std::sin(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */

		double i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
		double i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
		value *= - (i2 - i1) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double MissileField::electric_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

// Magnetic field

double MissileField::magnetic_rho (double vt, double rho, double phi, double z) const
{
	double value = this->A0 / 2;
	double R = this->R;
	mpf_class vt_z = mpf_class(vt * vt);
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		throw std::logic_error("MissileField::magnetic_rho is not implemented for -DNUMERIC_PDISK_LINEAR_INT");

	#else /* NUMERIC_PDISK_LINEAR_INT */
		
		value *= MissileField::int_lommel_001(vt, rho, z, R) - 
				 MissileField::int_lommel_011(vt, rho, z, R);
		return value * std::sin(phi);

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double MissileField::magnetic_phi (double vt, double rho, double phi, double z) const
{
	double value = this->A0 / 2;
	double R = this->R;
	mpf_class vt_z = mpf_class(vt * vt);
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		throw std::logic_error("MissileField::magnetic_phi is not implemented for -DNUMERIC_PDISK_LINEAR_INT");

	#else /* NUMERIC_PDISK_LINEAR_INT */

		value *= MissileField::int_lommel_011(vt, rho, z, R) * std::cos(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double MissileField::magnetic_z (double vt, double rho, double phi, double z) const
{
	double value = this->A0;
	double R = this->R;
	mpf_class vt_z = mpf_class(vt * vt);
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	if (vt_z.get_d() <= 0) return 0;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		throw std::logic_error("MissileField::magnetic_z is not implemented for -DNUMERIC_PDISK_LINEAR_INT");

	#else /* NUMERIC_PDISK_LINEAR_INT */

		value *= MissileField::int_lommel_111(vt, rho, z, R) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

// =============================================================================

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
	if (z + rho < 0.05) return this->A0/4;
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

double MissileField::int_lommel_001 (double vt, double rho, double z, double R)
{
	mp_bitcnt_t bitrate = mp_bitcnt_t(256);
	mpf_set_default_prec(bitrate);

	mpf_class ctpz(vt), ctmz(vt), ctz;
	mpf_add(ctpz.get_mpf_t(), ctpz.get_mpf_t(), mpf_class(z).get_mpf_t());
	mpf_add(ctmz.get_mpf_t(), ctmz.get_mpf_t(), mpf_class(-z).get_mpf_t());
	mpf_div(ctz.get_mpf_t(), ctmz.get_mpf_t(), ctpz.get_mpf_t());
	mpf_class ctz_in_pow(ctz);

	mpf_class ct_z = mpf_class(vt * vt);
	mpf_add( ct_z.get_mpf_t(), ct_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	
	mpf_class yacob_arg, yacob;

	if (rho == 0) {
		if (R < std::sqrt(ct_z.get_d())) {
			mpf_class sum = mpf_class(0);
			for (std::size_t m = 0; m <= MissileField::STATIC_TERMS_NUMBER; m++) {
				yacob = Math::yacobi_polinom(m, 1 - 2 * R * R / ct_z);
				mpf_mul(yacob.get_mpf_t(), yacob.get_mpf_t(), ctz_in_pow.get_mpf_t());
				mpf_add(sum.get_mpf_t(), sum.get_mpf_t(), yacob.get_mpf_t());
				mpf_mul(ctz_in_pow.get_mpf_t(), ctz_in_pow.get_mpf_t(), ctz.get_mpf_t());
			}
			mpf_class res = mpf_class(2 * R * R);
			mpf_mul(res.get_mpf_t(), res.get_mpf_t(), sum.get_mpf_t());
			mpf_div(res.get_mpf_t(), res.get_mpf_t(), ct_z.get_mpf_t());
			return res.get_d();
		}
		if (R >= std::sqrt(ct_z.get_d())) return 1;
	}

	throw std::logic_error("MissileField::int_lommel_001 is not implemented for non zero rho!");
}

double MissileField::int_lommel_011 (double vt, double rho, double z, double R)
{
	return MissileField::int_lommel_001(vt, rho, z, R) / 2;
}

double MissileField::int_lommel_111 (double vt, double rho, double z, double R)
{
	UNUSED(vt); UNUSED(rho); UNUSED(z); UNUSED(R);
	throw std::logic_error("MissileField::int_lommel_111 is not implemented");
}
