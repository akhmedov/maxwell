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

#include "maxwell.hpp"
#include "uniform_disk_current.hpp"
#include "updisk_meandr.hpp"

#include <regex>
#include <string>
#include <cmath>

struct KerrAmendment : public TransientResponse {

	KerrAmendment (double R, double A0, double eps_r, double mu_r, double chi3, Logger* global_logger = NULL);

	double electric_x (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double electric_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_rho (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_phi (const Point::SpaceTime<Point::Cylindrical>& event) const override;
	double magnetic_z (const Point::SpaceTime<Point::Cylindrical>& event) const override;

	double observed_from (const Point::Cylindrical& point) const override;
	double observed_to (const Point::Cylindrical& point) const override;

	double modal_jm  (const Point::ModalSpaceTime<Point::Cylindrical>& event) const override;
	double modal_vmh (const Point::ModalSpaceTime<Point::Cylindrical>& event) const override;

protected:

	static double alpha  (double ct_z, double rho, double R);
	static double beta   (double ct_z, double rho, double R);
	static double gamma  (double ct_z, double rho, double R);
	static double lambda (double ct_z, double rho, double R);

	static double int_bessel_011_perp (double ct_z, double rho, double R);
	static double int_bessel_001_perp (double ct_z, double rho, double R);

private:
	double kerr;
};

#endif /* kerr_amendment_hpp */
