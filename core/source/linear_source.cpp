//
//  linear_source.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 26.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "linear_source.hpp"

// Current distribution density

LinearSource::LinearSource (double tau)
: duration(tau) { }

double LinearSource::get_duration ()
{
	return this->duration;
}

void LinearSource::set_duration (double tau0)
{
	this->duration = tau0;
}

LinearCurrent::LinearCurrent (double tau)
: LinearSource(tau) { }

double LinearCurrent::x  (double rho, double phi, double z) const
{
	const double current_rho = this->rho(rho,phi,z);
	const double current_phi = this->phi(rho,phi,z);
	const double current_x = current_rho * std::cos(phi) - current_phi * std::sin(phi);
	return current_x;
}

double LinearCurrent::y  (double rho, double phi, double z) const
{
	const double current_rho = this->rho(rho,phi,z);
	const double current_phi = this->phi(rho,phi,z);
	const double current_y = current_rho * std::sin(phi) + current_phi * std::cos(phi);
	return current_y;
}

std::vector<double> LinearCurrent::cylindric (double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->rho(rho,phi,z));
	array.push_back(this->phi(rho,phi,z));
	array.push_back(this->z(rho,phi,z));
	return array;
}

std::vector<double> LinearCurrent::cartesian (double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->x(rho,phi,z));
	array.push_back(this->y(rho,phi,z));
	array.push_back(this->z(rho,phi,z));
	return array;
}

// Charge distribution density

LinearCharge::LinearCharge (double tau)
: LinearSource(tau) { }
