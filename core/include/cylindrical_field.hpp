//
//  cylindrical_field.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 24.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "integral.hpp"
#include "abstract_field.hpp"

template <class System> struct CylindricalField : public AbstractField<System> {

	using AbstractField<System>::AbstractField;

    virtual double electric_rho (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_phi (const Point::SpaceTime<System>& event) const = 0;
	virtual double electric_z (const Point::SpaceTime<System>& event) const = 0;

	virtual double magnetic_rho (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_phi (const Point::SpaceTime<System>& event) const = 0;
	virtual double magnetic_z (const Point::SpaceTime<System>& event) const = 0;

	virtual double observed_from (const System& point) const = 0;
	virtual double observed_to   (const System& point) const = 0;

    virtual double electric_x (const Point::SpaceTime<System>& event) const;
	virtual double electric_y (const Point::SpaceTime<System>& event) const;

	virtual double electric_r (const Point::SpaceTime<System>& event) const;
	virtual double electric_theta (const Point::SpaceTime<System>& event) const;

	virtual double magnetic_x (const Point::SpaceTime<System>& event) const;
	virtual double magnetic_y (const Point::SpaceTime<System>& event) const;

	virtual double magnetic_r (const Point::SpaceTime<System>& event) const;
	virtual double magnetic_theta (const Point::SpaceTime<System>& event) const;

	virtual double energy   (const System& point) const;
	virtual double energy_e (const System& point) const;
	virtual double energy_h (const System& point) const;
};

template <class System> double CylindricalField<System>::electric_x (const Point::SpaceTime<System>& event) const
{
    double phi = Point::Cylindrical::convert(event).phi();
	double Erho = this->electric_rho(event);
	double Ephi = this->electric_phi(event);
	return Erho * std::cos(phi) - Ephi * std::sin(phi);
}

template <class System> double CylindricalField<System>::electric_y (const Point::SpaceTime<System>& event) const
{
    double phi = Point::Cylindrical::convert(event).phi();
	double Erho = electric_rho(event);
	double Ephi = electric_phi(event);
	return Erho * std::sin(phi) + Ephi * std::cos(phi);
}

template <class System> double CylindricalField<System>::electric_r (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CylindricalField::electric_r is not implemented!");
}

template <class System> double CylindricalField<System>::electric_theta (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CylindricalField::electric_theta is not implemented!");
}

template <class System> double CylindricalField<System>::magnetic_x (const Point::SpaceTime<System>& event) const
{
    double phi = Point::Cylindrical::convert(event).phi();
	double Hrho = this->magnetic_rho(event);
	double Hphi = this->magnetic_phi(event);
	return Hrho * std::cos(phi) - Hphi * std::sin(phi);
}

template <class System> double CylindricalField<System>::magnetic_y (const Point::SpaceTime<System>& event) const
{
    double phi = Point::Cylindrical::convert(event).phi();
	double Hrho = magnetic_rho(event);
	double Hphi = magnetic_phi(event);
	return Hrho * std::sin(phi) + Hphi * std::cos(phi);
}

template <class System> double CylindricalField<System>::magnetic_r (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CylindricalField::magnetic_r is not implemented!");
}

template <class System> double CylindricalField<System>::magnetic_theta (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CylindricalField::magnetic_theta is not implemented!");
}

template <class System> double CylindricalField<System>::energy (const System&) const
{
    throw std::logic_error("CylindricalField::energy is not implemented!");
}

template <class System> double CylindricalField<System>::energy_e (const System& point) const
{
    double from = this->observed_from(point);
    double to = this->observed_to(point);

	auto f = [this, point] (double vt) {
		Point::SpaceTime<System> event{point};
        event.ct() = vt;
		double Erho = this->electric_rho(event);
		double Ephi = this->electric_phi(event);
		double Ez   = this->electric_z(event);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try {
		return SimpsonRunge(5e1, this->error, 1e5).value(from, to, f);
	} catch (double not_trusted) { // TODO: double?
		if (this->global_log) {
			std::string mesg = CylindricalField::INTEGRAL_WARNING;
			mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "CylindricalField::energy_e");
			mesg = std::regex_replace(mesg, std::regex("\\$POINT" ), point.to_str());
			this->global_log->warning(mesg);
		}
		return not_trusted;
	}
}

template <class System> double CylindricalField<System>::energy_h (const System& point) const
{
    double from = this->observed_from(point);
    double to = this->observed_to(point);

	auto f = [this, point] (double vt) {
		Point::SpaceTime<System> event{point};
        event.ct() = vt;
		double Erho = this->magnetic_rho(event);
		double Ephi = this->magnetic_phi(event);
		double Ez   = this->magnetic_z(event);
		return Erho*Erho + Ephi*Ephi + Ez*Ez;
	};

	try {
		return SimpsonRunge(5e1, this->error, 1e5).value(from, to, f);
	} catch (double not_trusted) { // TODO: double?
		if (this->global_log) {

			std::string mesg = CylindricalField::INTEGRAL_WARNING;
			mesg = std::regex_replace(mesg, std::regex("\\$NAME" ), "CylindricalField::energy_h");
			mesg = std::regex_replace(mesg, std::regex("\\$RHO" ), point.to_str());
			this->global_log->warning(mesg);
		}
		return not_trusted;
	}
}
