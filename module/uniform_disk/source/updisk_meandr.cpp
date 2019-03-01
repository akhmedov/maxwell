//
//  updisk_meandr.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "updisk_meandr.hpp"

// Linear current distribution of time squared shape

MeandrPeriod::MeandrPeriod (double disk_radius, double magnitude, double duration)
: UniformPlainDisk(disk_radius, magnitude) 
{
	this->tau = duration;
}

double MeandrPeriod::rho (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("MeandrPeriod::rho is not implemented");
}

double MeandrPeriod::phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("MeandrPeriod::phi is not implemented");
}

double MeandrPeriod::z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("MeandrPeriod::z is not implemented");
}

double MeandrPeriod::get_duration () const
{
	return this->tau;
}

UniformPlainDisk* MeandrPeriod::updisk () const
{
	return new UniformPlainDisk(this->R, this->A0);
}

// Electric field

SquaredPulse::SquaredPulse (MeandrPeriod* source, Homogeneous* medium)
: MissileField(source, medium) {
	this->A0 = source->get_magnitude();
	this->R = source->get_disk_radius();
	this->tau = source->get_duration();
}

double SquaredPulse::electric_rho (double vt, double rho, double phi, double z) const
{
	double R = this->R;
	double tau = this->tau;
	double value = this->A0 / 2;

	double sqrt_vt_z  = std::sqrt(vt*vt - z*z);
	double sqrt_tau_z = std::sqrt((vt-tau)*(vt-tau) - z*z);
	if (std::isnan(sqrt_vt_z) || (sqrt_vt_z == 0)) return 0;

	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(vt,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(vt,z));

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
			i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
		} else {
			i1 = SquaredPulse::int_bessel_011(sqrt_vt_z, sqrt_tau_z, rho, R);
		}

		return value * i1 * std::cos(phi);

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double SquaredPulse::electric_phi (double vt, double rho, double phi, double z) const
{
	double R = this->R;
	double tau = this->tau;
	double value = this->A0 / 2;

	double sqrt_vt_z  = std::sqrt(vt*vt - z*z);
	double sqrt_tau_z = std::sqrt((vt-tau)*(vt-tau) - z*z);
	if (std::isnan(sqrt_vt_z) || (sqrt_vt_z == 0)) return 0;

	value *= std::sqrt(LinearMedium::MU0 * medium->relative_permeability(vt,z));
	value /= std::sqrt(LinearMedium::EPS0 * medium->relative_permittivity(vt,z));

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
			i1 = MissileField::int_bessel_011(sqrt_vt_z, rho, R);
			i2 = MissileField::int_bessel_001(sqrt_vt_z, rho, R);
		} else {
			i1 = SquaredPulse::int_bessel_011(sqrt_vt_z, sqrt_tau_z, rho, R);
			i2 = SquaredPulse::int_bessel_001(sqrt_vt_z, sqrt_tau_z, rho, R);	
		}

		value *= - (i2 - i1) * std::sin(phi);
		return value;

	#endif /* NUMERIC_PDISK_LINEAR_INT */
}

double SquaredPulse::electric_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("SquaredPulse::electric_z is not implemented");
}

double SquaredPulse::magnetic_rho (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("SquaredPulse::magnetic_rho is not implemented");
}

double SquaredPulse::magnetic_phi (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("SquaredPulse::magnetic_phi is not implemented");
}

double SquaredPulse::magnetic_z (double vt, double rho, double phi, double z) const
{
	UNUSED(vt); UNUSED(rho); UNUSED(phi); UNUSED(z);
	throw std::logic_error("SquaredPulse::magnetic_z is not implemented");
}

double SquaredPulse::int_bessel_001 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R)
{
	double i1_from = MissileField::int_bessel_001(sqrt_vt_z,rho,R);
	double i1_to = MissileField::int_bessel_001(sqrt_tau_z,rho,R);
	return i1_from - i1_to;
}

double SquaredPulse::int_bessel_011 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R)
{
	double i1_from = MissileField::int_bessel_011(sqrt_vt_z,rho,R);
	double i1_to = MissileField::int_bessel_011(sqrt_tau_z,rho,R);
	return i1_from - i1_to;
}

// extern "C"  {

// 	const char* load_module_name ()
// 	{
// 		return "UniformDisk.MeandrMonocycle";
// 	}

// 	LinearMedium* load_module_medium ()
// 	{
// 		return new Homogeneous(1,1);
// 	}

// 	LinearCurrent* load_module_source ()
// 	{
// 		return new MeandrPeriod(1,1,0.5);
// 	}

// 	AbstractField* load_module_filed ()
// 	{
// 		auto medium = new Homogeneous(1,1);
// 		auto source = new MeandrPeriod(1,1,0.5);
// 		return new SquaredPulse(source,medium);
// 	}
// }
