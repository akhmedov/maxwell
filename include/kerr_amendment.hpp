//
//  kerr_amendment.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while ( false )
#endif

#ifndef kerr_amendment_hpp
#define kerr_amendment_hpp

#define NU_POINTS  2e3
#define VT_POINTS  20
#define PHO_POINTS 20
#define Z_POINTS   20

#include "integral.hpp"
#include "phys_math.hpp"
#include "nonlinear_medium.hpp"
#include "nonlinear_field.hpp"
#include "uniform_disk_current.hpp"

#include <complex>

using namespace std::complex_literals;

// TODO: rename class with name of medium that provides Kerr Effect
// TODO: implement Multiple inheritance (virtual inheritance)
struct KerrMedium : public NonlinearMedium, public Homogeneous {

	using Homogeneous::relative_permittivity;
	using Homogeneous::relative_permeability;
	
	KerrMedium (double realative_mu, double realative_eps, double xi3, double sigma);
	double conductivity (double ct, double z) const;
	double relative_permittivity (double ct, double z, std::size_t term) const;
	double relative_permeability (double ct, double z, std::size_t term) const;

private:
	/* double eps_r;
	double mu_r; */
	double sigma;
	double xi3;
};

struct KerrAmendment : public NonlinearField {

	KerrAmendment (MissileField* field, KerrMedium* medium, UniformPlainDisk* source);
	double electric_x (double ct, double rho, double phi, double z) const;

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

protected:

	static double x_trans (int m, double nu, double rho, double phi);

	static double N_sum (double R, int m, double nu, double ct, double rho_perp, double z); 
	static double N1    (double R, int m, double nu, double ct, double varrho, double z);
	static double N2    (double R, int m, double nu, double ct, double varrho, double z);
	static double N3    (double R, int m, double nu, double ct, double varrho, double z);
	static double N4    (double R, int m, double nu, double ct, double varrho, double z);

	static double int_bessel_011_perp (double vt, double z, double rho, double R);
	static double int_bessel_001_perp (double vt, double z, double rho, double R);

protected:
	
	MissileField* linear_field;
	double A0;
	double R;
};

#endif /* kerr_amendment_hpp */
