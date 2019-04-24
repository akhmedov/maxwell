//
//  updisk_meandr.cpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 13.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "updisk_meandr.hpp"

// Linear current distribution of time squared shape

SquaredPulse::SquaredPulse (double radius, double magnitude, double eps_r, double mu_r, double duration, Logger* global_log)
: TransientResponse(radius,magnitude,eps_r,mu_r,global_log), tau(duration) {}

double SquaredPulse::electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	Point::SpaceTime<Point::Cylindrical> event_tau = event; event_tau.ct() -= tau;
	double sqrt_tau_z = event_tau.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double rho = event.rho(), phi = event.phi();
	double value = A0 * std::sqrt(MU0 * MU / EPS0 * EPS) / 2;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		std::size_t bais = 1e5;
		Simpson I = Simpson(10*bais);

		auto ui_011 = [sqrt_vt_z, sqrt_tau_z, rho, R] (double nu) {
			if (nu == 0) return 0.0;
			double first, second;

			if (rho == 0) {
				first = R * j0(nu*sqrt_vt_z) * j1(nu*R) / 2;
				if (!std::isnan(sqrt_tau_z))
					second = R * j0(nu*sqrt_tau_z) * j1(nu*R) / 2;
				return first - second;
			}

			first = R * j0(nu*sqrt_vt_z) * j1(nu*rho) * j1(nu*R) / rho / nu;
			if (!std::isnan(sqrt_tau_z))
				second = R * j0(nu*sqrt_tau_z) * j1(nu*rho) * j1(nu*R) / rho / nu;
			return first - second;
		};

		value *= I.value(0, bais, ui_011) * std::cos(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */

		double i1;

		if (std::isnan(sqrt_tau_z) || (sqrt_tau_z == 0)) {
			i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
		} else {
			i1 = SquaredPulse::int_bessel_011(sqrt_vt_z, sqrt_tau_z, rho, R);
		}

		return value * i1 * std::cos(phi);

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double SquaredPulse::electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double sqrt_vt_z = event.sqrt_vt2_z2();
	Point::SpaceTime<Point::Cylindrical> event_tau = event; event_tau.ct() -= tau;
	double sqrt_tau_z = event_tau.sqrt_vt2_z2();
	if (std::isnan(sqrt_vt_z) || sqrt_vt_z == 0) return 0;

	double rho = event.rho(), phi = event.phi();
	double value = A0 * std::sqrt(MU0 * MU / EPS0 * EPS) / 2;

	#ifdef NUMERIC_PDISK_LINEAR_INT

		std::size_t bais = 1e5;
		Simpson I = Simpson(10*bais);

		auto ui_011 = [sqrt_vt_z, sqrt_tau_z, rho, R] (double nu) {
			if (nu == 0) return 0.0;
			double i1_first = 0, i1_second = 0;
			double i2_first = 0, i2_second = 0;

			if (rho == 0) {
				i1_first = j0(nu*sqrt_vt_z) * j1(nu*R) / 2;
				i2_first = j0(nu*sqrt_vt_z) * j1(nu*R);
				if (!std::isnan(sqrt_tau_z)) {
					i1_second = j0(nu*sqrt_tau_z) * j1(nu*R) / 2;
					i2_second = j0(nu*sqrt_tau_z) * j1(nu*R);
				}
				return R * (i1_first - i1_second - i2_first + i2_second);
			}

			i1_first = j0(nu*sqrt_vt_z) * j1(nu*rho) * j1(nu*R) / rho / nu;	
			i2_first = j0(nu*sqrt_vt_z) * j0(nu*rho) * j1(nu*R);	
			if (!std::isnan(sqrt_tau_z)) {
				i1_second = j0(nu*sqrt_tau_z) * j1(nu*rho) * j1(nu*R) / rho / nu;
				i2_second = j0(nu*sqrt_tau_z) * j0(nu*rho) * j1(nu*R);
			}

			return R * (i1_first - i1_second - i2_first + i2_second);
		};

		value *= I.value(0, bais, ui_011) * std::cos(phi);
		return value;

	#else /* NUMERIC_PDISK_LINEAR_INT */

		double i1 = 0, i2 = 0;
		if (std::isnan(sqrt_tau_z) || (sqrt_tau_z == 0)) {
			i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
			i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
		} else {
			i1 = SquaredPulse::int_bessel_011(sqrt_vt_z, sqrt_tau_z, rho, R);
			i2 = SquaredPulse::int_bessel_001(sqrt_vt_z, sqrt_tau_z, rho, R);	
		}

		value *= - (i2 - i1) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double SquaredPulse::electric_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	return 0;
}

double SquaredPulse::magnetic_rho (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("SquaredPulse::magnetic_rho is not implemented");
}

double SquaredPulse::magnetic_phi (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("SquaredPulse::magnetic_phi is not implemented");
}

double SquaredPulse::magnetic_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("SquaredPulse::magnetic_z is not implemented");
}

double SquaredPulse::int_bessel_001 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R)
{
	double i1_from = TransientResponse::int_bessel_001(sqrt_vt_z,rho,R);
	double i1_to = TransientResponse::int_bessel_001(sqrt_tau_z,rho,R);
	return i1_from - i1_to;
}

double SquaredPulse::int_bessel_011 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R)
{
	double i1_from = TransientResponse::int_bessel_011(sqrt_vt_z,rho,R);
	double i1_to = TransientResponse::int_bessel_011(sqrt_tau_z,rho,R);
	return i1_from - i1_to;
}

double SquaredPulse::observed_from (const Point::Cylindrical& point) const
{
	return TransientResponse::observed_from(point);
}

double SquaredPulse::observed_to (const Point::Cylindrical& point) const
{
	return this->tau + TransientResponse::observed_to(point);
}
