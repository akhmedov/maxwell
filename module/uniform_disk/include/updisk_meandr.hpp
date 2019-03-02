//
//  updisk_meandr.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#ifndef updisk_meandr_hpp
#define updisk_meandr_hpp

#include <gmp.h>
#include <gmpxx.h>

#include <exception>
#include <iostream>

#include "maxwell.hpp"
#include "uniform_disk_current.hpp"

struct MeandrPeriod : public UniformPlainDisk {

	MeandrPeriod (double disk_radius, double magnitude, double duration);
	double rho (double ct, double rho, double phi, double z) const;
	double phi (double ct, double rho, double phi, double z) const;
	double z (double ct, double rho, double phi, double z) const;

	// double get_magnitude () const;
	// double get_disk_radius () const;
	double get_duration () const;
	UniformPlainDisk* updisk () const;

protected:
	double A0;
	double R;
	double tau;
};

struct SquaredPulse : public MissileField {

	SquaredPulse (MeandrPeriod* source, Homogeneous* medium);

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

	static double int_bessel_001 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R); // I2 in thesis.pdf
	static double int_bessel_011 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R); // I1 in thesis.pdf
	
protected:

	double A0;
	double R;
	double tau;
};

#endif /* updisk_meandr_hpp */
