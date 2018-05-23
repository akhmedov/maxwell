//
//  electrodynamics.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "electrodynamics.hpp"

double Electrodynamics::directivity () const
{
	throw std::logic_error("Electrodynamics::directivity is not implemented!");
}

double Electrodynamics::energy_e (double rho, double phi, double z) const
{
	double min_vt = 0;
	double max_vt = 2 * std::sqrt(z*z + rho*rho);
	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	}; 
	return SimpsonRunge(1e3, 1, 1e5).value(min_vt,max_vt,f);
}

double Electrodynamics::energy_eh (double rho, double phi, double z) const
{
	double min_vt = 0;
	double max_vt = 2 * std::sqrt(z*z + rho*rho);
	auto f = [this, rho, phi, z] (double vt) {
		double Erho = this->electric_rho(vt,rho,phi,z);
		double Ephi = this->electric_phi(vt,rho,phi,z);
		double Ez = this->electric_z(vt,rho,phi,z);
		double Hrho = this->magnetic_rho(vt,rho,phi,z);
		double Hphi = this->magnetic_phi(vt,rho,phi,z);
		double Hz = this->magnetic_z(vt,rho,phi,z);
		return Erho*Hrho + Ephi*Hphi + Ez*Hz;
	}; 
	return SimpsonRunge(1e3, 1, 1e5).value(min_vt,max_vt,f);
}
