//
//  linear_duhamel.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 16.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "logger.hpp"
#include "integral.hpp"
#include "phys_math.hpp"
#include "abstract_field.hpp"

#include <regex>
#include <string>
#include <iostream>

#define INIT_NODES 1e3
#define MAX_NODES  1e6
#define STEP 	   1e-4

template <class System> using FieldComponent = std::function<double(AbstractField<System>* field, const Point::SpaceTime<System>& event)>;

template <template <class System> class AbstractFieldImpl, class System> struct DuhamelSuperpose : public AbstractFieldImpl<System> {

	DuhamelSuperpose (AbstractField<System>* tresponce, double duration, const std::function<double(double)>& shape, Logger* global_log = NULL)
	: AbstractFieldImpl<System>(global_log), tau0(duration), tr(tresponce), func(shape) { }

	double electric_x (const Point::SpaceTime<System>& event) const 
	{
		return this->duhamel(&AbstractField<System>::electric_x, event);
	}

	double electric_y (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_y, event);
	}

	double electric_z (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_z, event);
	}

	double electric_rho (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_rho, event);
	}

	double electric_phi (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_phi, event);
	}

	double electric_r (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_r, event);
	}

	double electric_theta (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::electric_theta, event);
	}

	double magnetic_x (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_x, event);
	}

	double magnetic_y (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_y, event);
	}

	double magnetic_z (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_z, event);
	}

	double magnetic_rho (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_rho, event);
	}

	double magnetic_phi (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_phi, event);
	}

	double magnetic_r (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_r, event);
	}

	double magnetic_theta (const Point::SpaceTime<System>& event) const
	{
		return this->duhamel(&AbstractField<System>::magnetic_theta, event);
	}

	double observed_from (const System& point) const
	{
		return this->tr->observed_from(point);
	}

	double observed_to (const System& point) const
	{
		return this->tau0 + this->tr->observed_to(point);
	}

	double duhamel (const FieldComponent<System>& comp, const Point::SpaceTime<System>& event) const
	{
		SimpsonRunge integral = SimpsonRunge(INIT_NODES, this->error, MAX_NODES);

		auto core = [this,&event,&comp] (double tau) {
			Point::SpaceTime<System> event_copy = event;
			event_copy.ct() -= tau;
			double shape = (tau == 0 || tau == this->tau0) ? 0 : Math::derivat4(func, tau);
			return shape * comp(this->tr,event_copy);
		};

		try {
			return integral.value(0, event.ct() - event.z(), core);
		} catch (double not_trusted) {

			if (this->global_log) {
				std::string mesg = DuhamelSuperpose::INTEGRAL_WARNING;
				mesg = std::regex_replace(mesg, std::regex("\\$NAME"), "DuhamelSuperpose<System>::duhamel");
				mesg = std::regex_replace(mesg, std::regex("\\$POINT"), event.to_str());
				this->global_log->warning(mesg);
			}

			Point::SpaceTime<System> prev{event};
			if (prev.ct() - STEP > 0) prev.ct() -= STEP;
			return comp(this->tr,prev);
		}
	}

private:

	double tau0;
	AbstractField<System>* tr;
	std::function<double(double)> func;
};
