//
//  abstract_field.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "abstract_field.hpp"

double AbstractField::electric_x (double ct, double rho, double phi, double z) const
{
	const double electric_rho = this->electric_rho(ct,rho,phi,z);
	const double electric_phi = this->electric_phi(ct,rho,phi,z);
	const double electric_x = electric_rho * std::cos(phi) - electric_phi * std::sin(phi);
	return electric_x;
}

double AbstractField::electric_y (double ct, double rho, double phi, double z) const
{
	const double electric_rho = this->electric_rho(ct,rho,phi,z);
	const double electric_phi = this->electric_phi(ct,rho,phi,z);
	const double electric_y = electric_rho * std::sin(phi) + electric_phi * std::cos(phi);
	return electric_y;
}

std::vector<double> AbstractField::electric_cylindric (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->electric_rho(ct, rho, phi, z));
	array.push_back(this->electric_phi(ct, rho, phi, z));
	array.push_back(this->electric_z(ct, rho, phi, z));
	return array;
}

std::vector<double> AbstractField::electric_cartesian (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->electric_x(ct, rho, phi, z));
	array.push_back(this->electric_y(ct, rho, phi, z));
	array.push_back(this->electric_z(ct, rho, phi, z));
	return array;
}
double AbstractField::magnetic_x (double ct, double rho, double phi, double z) const
{
	const double magnetic_rho = this->magnetic_rho(ct,rho,phi,z);
	const double magnetic_phi = this->magnetic_phi(ct,rho,phi,z);
	const double magnetic_x = magnetic_rho * std::cos(phi) - magnetic_phi * std::sin(phi);
	return magnetic_x;
}

double AbstractField::magnetic_y (double ct, double rho, double phi, double z) const
{
	const double magnetic_rho = this->magnetic_rho(ct,rho,phi,z);
	const double magnetic_phi = this->magnetic_phi(ct,rho,phi,z);
	const double magnetic_y = magnetic_rho * std::sin(phi) + magnetic_phi * std::cos(phi);
	return magnetic_y;
}

std::vector<double> AbstractField::magnetic_cylindric (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->magnetic_rho(ct, rho, phi, z));
	array.push_back(this->magnetic_phi(ct, rho, phi, z));
	array.push_back(this->magnetic_z(ct, rho, phi, z));
	return array;
}

std::vector<double> AbstractField::magnetic_cartesian (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.push_back(this->magnetic_x(ct, rho, phi, z));
	array.push_back(this->magnetic_y(ct, rho, phi, z));
	array.push_back(this->magnetic_z(ct, rho, phi, z));
	return array;
}
