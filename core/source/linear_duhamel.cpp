//
//  linear_duhamel.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 16.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "linear_duhamel.hpp"

const std::string LinearDuhamel::WARNING_MSG = "$COMP is not trusted at $TIME, $RHO, $PHI, $Z.";

FreeTimeCurrent::FreeTimeCurrent(LinearCurrent* on_source)
: base(on_source) { }

void FreeTimeCurrent::set_time_depth (const std::function<double(double)>& func)
{
	this->time_fnc = func;
}

double FreeTimeCurrent::time_shape (double vt) const
{
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

LinearDuhamel::LinearDuhamel (LinearCurrent* source, LinearMedium* medium, LinearField* on, Logger* log)
: LinearField(source,medium) 
{
	this->on_field = on;
	this->global_log = log;
}

void LinearDuhamel::set_accuracy (double persent)
{
	this->accuracy = persent;
}

double LinearDuhamel::electric_rho (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->electric_rho(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Erho");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->electric_rho (vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::electric_phi (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->electric_phi(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Ephi");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->electric_phi(vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::electric_z (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->electric_z(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Ez");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->electric_z(vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::magnetic_rho (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->magnetic_rho(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Hrho");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->magnetic_rho (vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::magnetic_phi (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->magnetic_phi(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Hphi");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->magnetic_phi (vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::magnetic_z (double vt, double rho, double phi, double z) const
{
	SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->accuracy, MAX_NODES);
	double tau0 = this->source->get_duration();

	auto core = [this,vt,tau0,rho,phi,z] (double tau) {
		auto f = [this,vt,rho,phi,z] (double tau) {
			return this->source->time_shape(tau);
		};
		double on_perp = this->on_field->magnetic_z(vt - tau,rho,phi,z);
		double shape = (tau == 0 || tau == tau0) ? 0 : Math::derivat4(f, tau);
		return shape * on_perp;
	};

	try { 
		return integral.value(0, vt-z, core);
	} catch (double not_trusted) {
		if (this->global_log) {
			std::string mesg = LinearDuhamel::WARNING_MSG;
			mesg = std::regex_replace(mesg, std::regex("\\$COMP"), "Hz");
			mesg = std::regex_replace(mesg, std::regex("\\$TIME"), std::to_string(vt));
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), std::to_string(rho));
			mesg = std::regex_replace(mesg, std::regex("\\$PHI" ), std::to_string(180*phi/M_PI));
			mesg = std::regex_replace(mesg, std::regex("\\$Z"   ), std::to_string(z));
			this->global_log->warning(mesg);
		}
		return this->magnetic_z (vt-STEP, rho, phi, z);
	}
}

double LinearDuhamel::observed_from (double rho, double phi, double z) const
{
	double real = this->on_field->observed_from(rho,phi,z);
	return (real - OBSERVAION_EPS > 0) ? real - OBSERVAION_EPS : 0;
}

double LinearDuhamel::observed_to (double rho, double phi, double z) const
{
	return this->source->get_duration() + this->on_field->observed_to(rho,phi,z) + OBSERVAION_EPS;
}
