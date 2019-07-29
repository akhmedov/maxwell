//
//  uniform_disk_current.cpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "uniform_disk_current.hpp"

// Electric field

std::size_t TransientResponse::STATIC_TERMS_NUMBER = 100;

TransientResponse::TransientResponse (double radius, double magnitude, double eps_r, double mu_r, Logger* global_log)
: CylindricalField<Point::Cylindrical>(global_log,1), A0(magnitude), R(radius), MU(mu_r), EPS(eps_r) {}

void TransientResponse::set_yterms_num (std::size_t number)
{
	TransientResponse::STATIC_TERMS_NUMBER = number;
}

std::size_t TransientResponse::get_yterms_num ()
{
	return TransientResponse::STATIC_TERMS_NUMBER;
}

double TransientResponse::electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double rho = event.rho(), phi = event.phi();
	double value = A0 * std::sqrt(MU0 * MU / EPS0 * EPS) / 2;

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

		double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
		value *= i1 * std::cos(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double TransientResponse::electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double rho = event.rho(), phi = event.phi();
	double value = A0 * std::sqrt(MU0 * MU / EPS0 * EPS) / 2;

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

		double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
		double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
		value *= - (i2 - i1) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double TransientResponse::electric_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	return 0;
}

// Magnetic field

double TransientResponse::magnetic_rho (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double value = A0 / 2;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		SimpsonRunge I = SimpsonRunge(10,1);

		auto ui_112 = [vt, rho, z, R, sqrt_vt_z] (double nu) {
			if (nu == 0) return 0.0;
			double W = - nu * (vt - z);
			double Z = nu * sqrt_vt_z;
			double U0 = Math::lommel (30, 2, W, Z); 
			if (rho == 0) return R * U0 * j1(nu*R) / 2;
			return j1(nu*rho) * j1(nu*R) * U0 / rho / nu;
		};

		double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
		double i3 = i1 + 2 * R * I.value(0, 1e4, ui_112);

		value *= i3 * std::cos(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */
		
		value *= TransientResponse::int_lommel_001(vt, rho, z, R) - 
				 TransientResponse::int_lommel_011(vt, rho, z, R);
		return value * std::sin(phi);

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double TransientResponse::magnetic_phi (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double value = A0 / 2;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		SimpsonRunge I = SimpsonRunge(10,1);

		auto ui_012 = [vt, rho, z, R, sqrt_vt_z] (double nu) {
			if (nu == 0) return 0.0;
			double W = - nu * (vt - z);
			double Z = nu * sqrt_vt_z;
			double U0 = Math::lommel (30, 2, W, Z); 
			return j0(nu*rho) * j1(nu*R) * U0; 
		};
		
		auto ui_112 = [vt, rho, z, R, sqrt_vt_z] (double nu) {
			if (nu == 0) return 0.0;
			double W = - nu * (vt - z);
			double Z = nu * sqrt_vt_z;
			double U0 = Math::lommel (30, 2, W, Z); 
			if (rho == 0) return R * U0 * j1(nu*R) / 2;
			return j1(nu*rho) * j1(nu*R) * U0 / rho / nu;
		};

		double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
		double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
		double i3 = i1 + 2 * R * I.value(0, 1e4, ui_112);
		double i4 = i2 + R * I.value(0, 1e4, ui_012);

		value *= - (i4 - i3) * std::sin(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */

		value *= TransientResponse::int_lommel_011(vt, rho, z, R) * std::cos(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double TransientResponse::magnetic_z (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double vt = event.ct(), rho = event.rho(), phi = event.phi(), z = event.z();
	double value = A0 / 2;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		throw std::logic_error("TransientResponse::magnetic_z is not implemented for -DNUMERIC_PDISK_LINEAR_INT");

	#else /* NUMERIC_PDISK_LINEAR_INT */

		value *= TransientResponse::int_lommel_111(vt, rho, z, R) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

// =============================================================================

double TransientResponse::static_magnitude (double z) const
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

double TransientResponse::static_magnitude (const Point::Cylindrical& point, double eps) const
{
	if (point.z() + point.rho() < 0.05) return this->A0/4;
	double R2 = this->R * this->R;
	double z2 = point.z()*point.z();
	double ct = point.rho() + std::sqrt(R2 + z2) + eps;
	Point::SpaceTime<Point::Cylindrical> event{point};
	event.ct() = ct;
	return this->magnetic_y(event);
}

// integrals

double TransientResponse::int_bessel_001 (double sqrt_vt_z, double rho, double R)
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

double TransientResponse::int_bessel_011 (double sqrt_vt_z, double rho, double R)
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

double TransientResponse::int_lommel_001 (double vt, double rho, double z, double R)
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
			for (std::size_t m = 0; m <= TransientResponse::STATIC_TERMS_NUMBER; m++) {
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

	throw std::logic_error("TransientResponse::int_lommel_001 is not implemented for non zero rho!");
}

double TransientResponse::int_lommel_011 (double vt, double rho, double z, double R)
{
	return TransientResponse::int_lommel_001(vt, rho, z, R) / 2;
}

double TransientResponse::int_lommel_111 (double /*vt*/, double /*rho*/, double /*z*/, double /*R*/)
{
	throw std::logic_error("TransientResponse::int_lommel_111 is not implemented");
}

double TransientResponse::observed_from (const Point::Cylindrical& point) const
{
	double rho = point.rho(), z = point.z();
	double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z)-0.01 : z-0.01;
	return (from > 0) ? from : 0.005;
}

double TransientResponse::observed_to (const Point::Cylindrical& point) const
{
	double rho = point.rho(), z = point.z();
	return std::sqrt((rho+R)*(rho+R) + z*z) + 0.01;
}
