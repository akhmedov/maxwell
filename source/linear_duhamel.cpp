//
//  linear_duhamel.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 16.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "linear_duhamel.hpp"

FreeTimeCurrent::FreeTimeCurrent(LinearCurrent* on_source, double tau)
: LinearCurrent(tau), base(on_source) { }

void FreeTimeCurrent::set_time_depth (const std::function<double(double)>& func)
{
	this->time_fnc = func;
}

double FreeTimeCurrent::time_shape (double vt) const
{
	if (vt < 0 || vt > this->duration) return 0;
	return this->time_fnc(vt);
}

double FreeTimeCurrent::rho (double rho, double phi, double z) const
{
	return this->base->rho(rho, phi, z);
}

double FreeTimeCurrent::phi (double rho, double phi, double z) const
{
	return this->base->phi(rho, phi, z);
}

double FreeTimeCurrent::z (double rho, double phi, double z) const
{
	return this->base->z(rho, phi, z);
}

// Linear Field by Duramel integral

LinearDuramel::LinearDuramel (LinearCurrent* source, LinearMedium* medium, LinearField* on, Logger* log)
: LinearField(source,medium) 
{
	this->on_field = on;
	this->global_log = log;
}

void LinearDuramel::set_accuracy (double persent)
{
	this->accuracy = persent;
}

double LinearDuramel::electric_rho (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,rho,phi,z] (double tau) {
		auto field = [this,vt,rho,phi,z] (double tau) {
			return this->on_field->electric_rho(vt - tau,rho,phi,z); 
		};
		double on_perp = (tau != 0) ? Math::derivat4(field, tau) : 0;
		double res = this->source->time_shape(tau) * on_perp;
		return res;
	};

	try { 
		return integral.value(0, tau0, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string msg = "Erho is not trusted at " + std::to_string(vt);
			this->global_log->warning(msg);
		}
		return not_trusted;
	}
}

double LinearDuramel::electric_phi (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,rho,phi,z] (double tau) {
		auto field = [this,vt,rho,phi,z] (double tau) {
			return this->on_field->electric_phi(vt - tau,rho,phi,z); 
		};	
		double on_perp = (tau != 0) ? Math::derivat4(field, tau) : 0;
		double res = this->source->time_shape(tau) * on_perp;
		return res;
	};

	try { 
		return integral.value(0, tau0, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string msg = "Ephi is not trusted at " + std::to_string(vt);
			this->global_log->warning(msg);
		}
		return not_trusted;
	}
}

double LinearDuramel::electric_z (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("LinearDuramel::electric_z is not implemented");
}

double LinearDuramel::magnetic_rho (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("LinearDuramel::electric_z is not implemented");
}

double LinearDuramel::magnetic_phi (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("LinearDuramel::electric_z is not implemented");
}

double LinearDuramel::magnetic_z (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("LinearDuramel::electric_z is not implemented");
}
