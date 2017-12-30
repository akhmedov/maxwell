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

#include "config.hpp"
#include "integral.hpp"
#include "phys_math.hpp"
#include "linear_medium.hpp"
#include "linear_source.hpp"
#include "linear_field.hpp"

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
	double rho (double ct, double rho, double phi, double z) const;
	double phi (double ct, double rho, double phi, double z) const;
	double z (double ct, double rho, double phi, double z) const;

	double get_magnitude () const;
	double get_disk_radius () const;

protected:
	double A0;
	double R;
};

struct MissileField : public LinearField {

	MissileField (UniformPlainDisk* source, Homogeneous* medium);

	void set_yterms_num (std::size_t);
	std::size_t get_yterms_num () const;

	double static_magnitude (double z) const;
	double static_magnitude (double rho, double phi, double z, double eps = 10e-10) const;

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

protected:
	static double int_bessel_001 (double vt_z, double rho, double R); // I2 in thesis.pdf
	static double int_bessel_011 (double vt_z, double rho, double R); // I1 in thesis.pdf
	
	// TODO: refactor for integrtion test compatability
	double int_lommel_001 (double ct, double rho, double z) const; // I4 in thesis.pdf
	double int_lommel_011 (double ct, double rho, double z) const; // I3 in thesis.pdf
	double int_lommel_111 (double ct, double rho, double z) const; // I5 in thesis.pdf
	
	double A0;
	double R;
	std::size_t STATIC_TERMS_NUMBER = 100;
};

#endif /* uniform_disk_current_hpp */
