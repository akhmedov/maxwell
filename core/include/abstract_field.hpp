//
//  abstract_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include <cmath>
#include <string>
#include <vector>

// Forward declaration.
struct Logger;

struct AbstractField {
	AbstractField(Logger *log = nullptr, double acc = 1)
		: global_log(log),
		  accuracy(acc)
        {}
	virtual ~AbstractField () = default;
	
	virtual double electric_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_z (double ct, double rho, double phi, double z) const = 0;

	virtual double electric_x (double ct, double rho, double phi, double z) const;
	virtual double electric_y (double ct, double rho, double phi, double z) const;
	std::vector<double> electric_cylindric (double ct, double rho, double phi, double z) const;
	std::vector<double> electric_cartesian (double ct, double rho, double phi, double z) const;

	virtual double magnetic_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_z (double ct, double rho, double phi, double z) const = 0;

	virtual double magnetic_x (double ct, double rho, double phi, double z) const;
	virtual double magnetic_y (double ct, double rho, double phi, double z) const;
	std::vector<double> magnetic_cylindric (double ct, double rho, double phi, double z) const;
	std::vector<double> magnetic_cartesian (double ct, double rho, double phi, double z) const;

	virtual double energy_cart (double x, double y, double z, double from, double to) const;
	virtual double energy (double rho, double phi, double z, double from, double to) const;

	virtual double observed_from (double x, double y, double z) const = 0;
	virtual double observed_to (double x, double y, double z) const = 0;

	constexpr static const double C = 299792458;
	constexpr static const double C2 = C * C;
	constexpr static const double EPS0 = 10e7 / (4 * M_PI * C2);
	constexpr static const double MU0 = 4 * M_PI * 10e-7;

protected:
	Logger *global_log{};
	double accuracy{};
};

struct ZeroField : public AbstractField {
	double electric_rho(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double electric_phi(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double electric_z(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double magnetic_rho(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double magnetic_phi(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double magnetic_z(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const override
        {
		return 0;
	}
	double energy_cart(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const
        {
		return 0;
	}
	double energy(double/*ct*/, double/*rho*/, double/*phi*/, double/*z*/) const 
        {
		return 0;
	}
};
