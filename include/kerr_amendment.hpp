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

#define TERMS_NUMBER 10e3

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
	
	KerrMedium (double realative_mu, double realative_eps, double chi3_electric, double sigma);
	double conductivity (double ct, double z) const;
	double relative_permittivity (double ct, double z, std::size_t term) const;
	double relative_permeability (double ct, double z, std::size_t term) const;

private:
	/* double eps_r;
	double mu_r; */
	double sigma;
	double chi3;
};

struct KerrAmendment : public MissileField, public NonlinearField {

	KerrAmendment (MissileField* field, KerrMedium* medium);

	std::complex<double> complex_electric_rho (double ct, double rho, double phi, double z) const;
	std::complex<double> complex_electric_phi (double ct, double rho, double phi, double z) const;
	std::complex<double> complex_electric_z (double ct, double rho, double phi, double z) const;

	std::complex<double> complex_magnetic_rho (double ct, double rho, double phi, double z) const;
	std::complex<double> complex_magnetic_phi (double ct, double rho, double phi, double z) const;
	std::complex<double> complex_magnetic_z (double ct, double rho, double phi, double z) const;

	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;

	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

protected:

	double im_evolution_h (long m, double nu, double vt, double z) const;
	double im_evolution_Ih (long m, double nu, double vt, double z) const;
	double im_evolution_Vh (long m, double nu, double vt, double z) const;
	double im_modal_source (long m, double nu, double vt, double z) const;
	
	double riemann (double nu, double vt_diff, double z_diff) const;

	double N1 (long m, double nu, double ct, double z) const;
	double N2 (long m, double nu, double ct, double z) const;
	double N3 (long m, double nu, double ct, double z) const;
	double N4 (long m, double nu, double ct, double z) const;
	double N5 (long m, double nu, double ct, double z) const;
	double N6 (long m, double nu, double ct,  double z) const;

	double vint_bessel_011_perp (double vt, double z, double rho, double R) const;
	double vint_bessel_001_perp (double vt, double z, double rho, double R) const;
};

#endif /* kerr_amendment_hpp */
