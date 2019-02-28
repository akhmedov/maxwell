//
//  abstract_field.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "abstract_field.hpp"

#include <cmath>
#include <regex>

#include "integral.hpp"
#include "logger.hpp"

namespace {
const char *exept_mgs = "Integral in $NAME is not trusted at $RHO, $PHI, $Z";
}

double AbstractField::electric_x (double ct, double rho, double phi, double z) const
{
	auto e_rho = electric_rho(ct, rho, phi, z);
	auto e_phi = electric_phi(ct, rho, phi, z);
	return e_rho * std::cos(phi) - e_phi * std::sin(phi);
}

double AbstractField::electric_y (double ct, double rho, double phi, double z) const
{
	auto e_rho = electric_rho(ct, rho, phi, z);
	auto e_phi = electric_phi(ct, rho, phi, z);
	return e_rho * std::sin(phi) + e_phi * std::cos(phi);
}

std::vector<double> AbstractField::electric_cylindric (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.reserve(3);

	array.push_back(electric_rho(ct, rho, phi, z));
	array.push_back(electric_phi(ct, rho, phi, z));
	array.push_back(electric_z(ct, rho, phi, z));

	return array;
}

std::vector<double> AbstractField::electric_cartesian (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.reserve(3);

	array.push_back(electric_x(ct, rho, phi, z));
	array.push_back(electric_y(ct, rho, phi, z));
	array.push_back(electric_z(ct, rho, phi, z));

	return array;
}

double AbstractField::magnetic_x (double ct, double rho, double phi, double z) const
{
	auto m_rho = magnetic_rho(ct,rho,phi,z);
	auto m_phi = magnetic_phi(ct,rho,phi,z);
	return m_rho * std::cos(phi) - m_phi * std::sin(phi);
}

double AbstractField::magnetic_y(double ct, double rho, double phi, double z) const {
	auto m_rho = magnetic_rho(ct,rho,phi,z);
	auto m_phi = magnetic_phi(ct,rho,phi,z);
	return m_rho * std::sin(phi) + m_phi * std::cos(phi);
}

std::vector<double> AbstractField::magnetic_cylindric (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.reserve(3);

	array.push_back(magnetic_rho(ct, rho, phi, z));
	array.push_back(magnetic_phi(ct, rho, phi, z));
	array.push_back(magnetic_z(ct, rho, phi, z));

	return array;
}

std::vector<double> AbstractField::magnetic_cartesian (double ct, double rho, double phi, double z) const
{
	std::vector<double> array;
	array.reserve(3);

	array.push_back(magnetic_x(ct, rho, phi, z));
	array.push_back(magnetic_y(ct, rho, phi, z));
	array.push_back(magnetic_z(ct, rho, phi, z));

	return array;
}

double AbstractField::energy (double rho, double phi, double z, double from, double to) const
{
	auto f = [this, rho, phi, z] (double vt) {
		double Erho = electric_rho(vt,rho,phi,z);
		double Ephi = electric_phi(vt,rho,phi,z);
		double Ez = electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try {
		return SimpsonRunge(5e1, accuracy, 1e5).value(from, to, f);
	} catch (double not_trusted) { // TODO: double?
		if (global_log) {
			auto mesg = std::regex_replace(exept_mgs, std::regex("\\$NAME" ), "AbstractField::energy");
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			global_log->warning(mesg);
		}
		return not_trusted;
	}
}

double AbstractField::energy_cart (double x, double y, double z, double from, double to) const
{
	double rho = std::sqrt(x*x + y*y);
	double phi = std::atan2(y, x);

	auto f = [this, rho, phi, z] (double vt) {
		double Erho = electric_rho(vt,rho,phi,z);
		double Ephi = electric_phi(vt,rho,phi,z);
		double Ez = electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try {
		return SimpsonRunge(5e1, accuracy, 1e5).value(from, to,f);
	} catch (double not_trusted) { // TODO: double?
		if (global_log) {
			auto mesg = std::regex_replace(exept_mgs, std::regex("\\$NAME" ), "AbstractField::energy_cart");
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			global_log->warning(mesg);
		}
		return not_trusted;
	}
}
