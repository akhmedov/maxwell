//
//  uniform_disk_current.hpp
//  Evolution
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#ifndef uniform_disk_current_hpp
#define uniform_disk_current_hpp

#include <gmp.h>
#include <gmpxx.h>

#include <exception>
#include <iostream>

#include "integral.hpp"
#include "phys_math.hpp"
#include "linear_field.hpp"
#include "linear_medium.hpp"
#include "linear_source.hpp"

struct Homogeneous : public LinearMedium {

	using LinearMedium::relative_permittivity;
	using LinearMedium::relative_permeability;
	
	Homogeneous(double realative_mu, double realative_eps);
	double relative_permittivity (double ct, double z) const;
	double relative_permeability (double ct, double z) const;

protected:
	double mu_r;
	double eps_r;	
};

struct UniformPlainDisk : public LinearCurrent {

	UniformPlainDisk (double disk_radius, double magnitude);
	double time_shape (double vt) const;
	double rho (double rho, double phi, double z) const;
	double phi (double rho, double phi, double z) const;
	double z (double rho, double phi, double z) const;

	double get_magnitude () const;
	double get_disk_radius () const;

protected:
	double A0;
	double R;
};

struct MissileField : public LinearField {

	MissileField (UniformPlainDisk* source, Homogeneous* medium);

	static void set_yterms_num (std::size_t);
	static std::size_t get_yterms_num ();

	double static_magnitude (double z) const;
	double static_magnitude (double rho, double phi, double z, double eps = 10e-10) const;

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

	// double energy_cart (double x, double y, double z) const;

	static double int_bessel_001 (double vt_z, double rho, double R); // I2 in thesis.pdf
	static double int_bessel_011 (double vt_z, double rho, double R); // I1 in thesis.pdf

	static double int_lommel_001 (double ct, double rho, double z, double R); // I4 in thesis.pdf
	static double int_lommel_011 (double ct, double rho, double z, double R); // I3 in thesis.pdf
	static double int_lommel_111 (double ct, double rho, double z, double R); // I5 in thesis.pdf
	
protected:

	double A0;
	double R;
	static std::size_t STATIC_TERMS_NUMBER;
};

#endif /* uniform_disk_current_hpp */
