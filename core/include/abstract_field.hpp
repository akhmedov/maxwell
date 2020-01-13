//
//  abstract_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "space_point.hpp"

#include <cmath>
#include <string>
#include <vector>
#include <limits>
#include <regex>

// Forward declaration.
struct Logger;

template <class System> struct AbstractField {
	
    AbstractField (Logger* log, double err = 3 /* % */) : error(err), global_log(log) {}
	virtual ~AbstractField () = default;

	virtual double electric_x (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_y (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_z (const Point::SpaceTime<System>& event) const = 0;

	virtual double electric_rho (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_phi (const Point::SpaceTime<System>& event) const = 0;

	virtual double electric_r (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_theta (const Point::SpaceTime<System>& event) const = 0;

	virtual double magnetic_x (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_y (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_z (const Point::SpaceTime<System>& event) const = 0;

	virtual double magnetic_rho (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_phi (const Point::SpaceTime<System>& event) const = 0;

	virtual double magnetic_r (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_theta (const Point::SpaceTime<System>& event) const = 0;

	virtual double energy   (const System& point) const = 0;
	virtual double energy_e (const System& point) const = 0;
	virtual double energy_h (const System& point) const = 0;

	virtual double observed_from (const System& point) const = 0;
	virtual double observed_to   (const System& point) const = 0;

	constexpr static const double C = 299792458;
	constexpr static const double C2 = C * C;
	constexpr static const double EPS0 = 10e7 / (4 * M_PI * C2);
	constexpr static const double MU0 = 4 * M_PI * 10e-7;

protected:

    double error; // %
    Logger* global_log{};
    static const std::string INTEGRAL_WARNING;
};

template <class System> const std::string AbstractField<System>::INTEGRAL_WARNING = "Integral in $NAME is not trusted at $POINT";

template <class System> struct ZeroField : public AbstractField<System> {
	ZeroField() : AbstractField<System>(NULL,0) {}
	double electric_x (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_y (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_z (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_rho (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_phi (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_r (const Point::SpaceTime<System>&) const override { return 0; }
	double electric_theta (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_x (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_y (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_z (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_rho (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_phi (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_r (const Point::SpaceTime<System>&) const override { return 0; }
	double magnetic_theta (const Point::SpaceTime<System>&) const override { return 0; }
	double energy   (const System&) const override { return 0; }
	double energy_e (const System&) const override { return 0; }
	double energy_h (const System&) const override { return 0; }
	double observed_from (const System&) const override { return 0; }
	double observed_to   (const System&) const override { return 0; }
};
