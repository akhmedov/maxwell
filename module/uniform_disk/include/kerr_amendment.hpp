//
//  kerr_amendment.hpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 22.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while ( false )
#endif

#ifndef kerr_amendment_hpp
#define kerr_amendment_hpp

#define PERIODS_NU	3
#define NODES_NU	7
#define MIN_NODES	5
#define MAX_ERROR	10 // %

#include "maxwell.hpp"
#include "uniform_disk_current.hpp"

#include <regex>
#include <string>
#include <complex>

using namespace std::complex_literals;

struct KerrAmendment : protected TransientResponse {

	KerrAmendment (double R, double A0, double eps_r, double mu_r, double chi3, Logger* global_logger = NULL);

	// TODO: move from here
	double current_x (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double current_rho (const Point::SpaceTime<Point::Cylindrical>& event) const;
	double current_phi (const Point::SpaceTime<Point::Cylindrical>& event) const;

	double electric_x (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double observed_from (const Point::Cylindrical& point) const override;
	double observed_to (const Point::Cylindrical& point) const override;

protected:

	static double x_trans (int m, double nu, double rho, double phi);

	static double N_sum (double R, int m, double nu, double ct, double varrho, double z); 
	static double N1    (double R, int m, double nu, double ct, double varrho, double z);
	static double N2    (double R, int m, double nu, double ct, double varrho, double z);
	static double N3    (double R, int m, double nu, double ct, double varrho, double z);
	static double N4    (double R, int m, double nu, double ct, double varrho, double z);

	static double int_bessel_011_perp (double vt, double z, double rho, double R);
	static double int_bessel_001_perp (double vt, double z, double rho, double R);

private:
	double kerr;
};

#endif /* kerr_amendment_hpp */
