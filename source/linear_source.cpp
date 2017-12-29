//
//  linear_source.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 26.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "linear_source.hpp"

double LinearCurrent::x  (double ct, double rho, double phi, double z) const
{
	const double current_rho = this->rho(ct,rho,phi,z);
	const double current_phi = this->phi(ct,rho,phi,z);
	const double current_x = current_rho * std::cos(phi) - current_phi * std::sin(phi);
	return current_x;
}

double LinearCurrent::y  (double ct, double rho, double phi, double z) const
{
	const double current_rho = this->rho(ct,rho,phi,z);
	const double current_phi = this->phi(ct,rho,phi,z);
	const double current_y = current_rho * std::sin(phi) + current_phi * std::cos(phi);
	return current_y;
}

std::vector<double> LinearCurrent::cylindric (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->rho(ct,rho,phi,z));
	array.push_back(this->phi(ct,rho,phi,z));
	array.push_back(this->z(ct,rho,phi,z));
	return array;
}

std::vector<double> LinearCurrent::cartesian (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->x(ct,rho,phi,z));
	array.push_back(this->y(ct,rho,phi,z));
	array.push_back(this->z(ct,rho,phi,z));
	return array;
}
