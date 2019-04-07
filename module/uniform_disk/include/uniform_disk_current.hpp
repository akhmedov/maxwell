//
//  uniform_disk_current.hpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef uniform_disk_current_hpp
#define uniform_disk_current_hpp

#include <gmp.h>
#include <gmpxx.h>

#include <exception>
#include <iostream>

#include "maxwell.hpp"

struct TransientResponse : public CylindricalField<Point::Cylindrical> {

	TransientResponse (double R, double A0, double eps_r, double mu_r, Logger* global_log);

	static void set_yterms_num (std::size_t);
	static std::size_t get_yterms_num ();

	double static_magnitude (double z) const;
	double static_magnitude (const Point::Cylindrical& event, double eps = 10e-10) const;

	double electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;

	double magnetic_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;

	double observed_from (const Point::Cylindrical& point) const override;
	double observed_to (const Point::Cylindrical& point) const override;

	static double int_bessel_001 (double vt_z, double rho, double R); // I2 in thesis.pdf
	static double int_bessel_011 (double vt_z, double rho, double R); // I1 in thesis.pdf

	static double int_lommel_001 (double ct, double rho, double z, double R); // I4 in thesis.pdf
	static double int_lommel_011 (double ct, double rho, double z, double R); // I3 in thesis.pdf
	static double int_lommel_111 (double ct, double rho, double z, double R); // I5 in thesis.pdf

protected:

	double A0;
	double R;
	double MU;
	double EPS;
	static std::size_t STATIC_TERMS_NUMBER;
};

#endif /* uniform_disk_current_hpp */
