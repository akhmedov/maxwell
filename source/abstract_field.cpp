//
//  abstract_field.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "abstract_field.hpp"

const std::string AbstractField::int_exept_mgs = "Integral in $NAME is not trusted at $RHO, $PHI, $Z";

AbstractField::AbstractField (Logger* log, double acc)
: accuracy(acc), global_log(log) { }

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

double AbstractField::energy (double rho, double phi, double z) const
{
	double min_vt = 0;
	double tau1 = std::sqrt((rho-1)*(rho-1) + z*z);
	double tau2 = tau1 + 0.5;

	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try { 
		return SimpsonRunge(5e1, this->accuracy, 1e5).value(tau1,tau2,f); 
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = AbstractField::int_exept_mgs;
			mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "AbstractField::energy");
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return not_trusted; 
	}
}

double AbstractField::energy_cart (double x, double y, double z) const
{
	double R = 1;
	double tau = 0.5;

	double rho = std::sqrt(x*x + y*y);
	double phi = std::atan2(y, x);
	double tau1 = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
	double tau2 = tau1 + 2*tau;

	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try { 
		return SimpsonRunge(5e1, this->accuracy, 1e5).value(tau1,tau2,f); 
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = AbstractField::int_exept_mgs;
			mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "AbstractField::energy_cart");
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return not_trusted; 
	}
}
