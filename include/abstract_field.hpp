//
//  abstract_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef abstract_field_hpp
#define abstract_field_hpp

#include "logger.hpp"
#include "integral.hpp"

#include <cmath>
#include <string>
#include <vector>

struct AbstractField {

	AbstractField (Logger* global_log = NULL, double accuracy = 1);
	
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

	virtual double energy_cart (double x, double y, double z) const;
	virtual double energy (double rho, double phi, double z) const;

protected:

	static const std::string int_exept_mgs;
	Logger* global_log;
	double accuracy;
};

#endif /* abstract_field_hpp */
