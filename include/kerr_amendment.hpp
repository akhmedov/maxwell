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

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

protected:

	double im_modal_source_sum (double nu, double ct, double varrho, double z) const;
	double im_modal_source (int m, double nu, double ct, double varrho, double z) const;
	double riemann (double nu, double vt_diff, double z_diff) const; 

	double N1 (int m, double nu, double ct, double varrho, double z) const;
	double N2 (int m, double nu, double ct, double varrho, double z) const;
	double N4 (int m, double nu, double ct, double varrho, double z) const;
	double N5 (int m, double nu, double ct, double varrho, double z) const;

	double int_bessel_011_perp (double vt, double z, double rho, double R) const;
	double int_bessel_001_perp (double vt, double z, double rho, double R) const;

protected:
	
	MissileField* linear_field;
	double A0;
	double R;
};

#endif /* kerr_amendment_hpp */
