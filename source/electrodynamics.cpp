//
//  electrodynamics.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "electrodynamics.hpp"

double Electrodynamics::energy (double rho, double phi, double z) const
{
	double min_vt = 0;
	double max_vt = 2 * std::sqrt(z*z + rho*rho);
	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	}; 
	return SimpsonRunge(1e2, 1, 1e5).value(min_vt,max_vt,f);
}

double Electrodynamics::energy_cart (double x, double y, double z) const
{
	double rho = std::sqrt(x*x + y*y);
	double phi = std::atan2(y, x);
	double min_vt = 0;
	double max_vt = 2 * std::sqrt(z*z + rho*rho);
	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	}; 
	return SimpsonRunge(1e2, 1, 1e5).value(min_vt,max_vt,f);	
}