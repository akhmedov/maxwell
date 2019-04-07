//
//  updisk_meandr.hpp
//  uniform_disk.module.maxwell
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

struct SquaredPulse : public TransientResponse {

	SquaredPulse (double radius, double magnitude, double eps_r, double mu_r, double duration, Logger* global_log);

	double electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double electric_z (const Point::SpaceTime<Point::Cylindrical>& event) const;

	double magnetic_rho (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double magnetic_phi (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double magnetic_z (const Point::SpaceTime<Point::Cylindrical>& event) const;

	static double int_bessel_001 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R); // I2 in thesis.pdf
	static double int_bessel_011 (double sqrt_vt_z, double sqrt_tau_z, double rho, double R); // I1 in thesis.pdf
	
protected:

	double tau;
};

#endif /* updisk_meandr_hpp */
