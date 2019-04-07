//
//  cartesian_field.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 24.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "integral.hpp"
#include "abstract_field.hpp"

template <class System> struct CartesianField : public AbstractField<System> {

	using AbstractField<System>::AbstractField;

	virtual double electric_x (const Point::SpaceTime<System>& event) const = 0; // override? using?
	virtual double electric_y (const Point::SpaceTime<System>& event) const = 0; // override? using?
	virtual double electric_z (const Point::SpaceTime<System>& event) const = 0; // override? using?

	virtual double magnetic_x (const Point::SpaceTime<System>& event) const = 0; // override? using?
	virtual double magnetic_y (const Point::SpaceTime<System>& event) const = 0; // override? using?
	virtual double magnetic_z (const Point::SpaceTime<System>& event) const = 0; // override? using?

	virtual double observed_from (const System& point) const = 0; // override? using?
	virtual double observed_to   (const System& point) const = 0; // override? using?

	virtual double electric_rho (const Point::SpaceTime<System>& event) const override;
	virtual double electric_phi (const Point::SpaceTime<System>& event) const override;

	virtual double electric_r (const Point::SpaceTime<System>& event) const override;
	virtual double electric_theta (const Point::SpaceTime<System>& event) const override;

	virtual double magnetic_rho (const Point::SpaceTime<System>& event) const override;
	virtual double magnetic_phi (const Point::SpaceTime<System>& event) const override;

	virtual double magnetic_r (const Point::SpaceTime<System>& event) const override;
	virtual double magnetic_theta (const Point::SpaceTime<System>& event) const override;

	virtual double energy   (const System& point) const;
	virtual double energy_e (const System& point) const;
	virtual double energy_h (const System& point) const;
};

template <class System> double CartesianField<System>::electric_rho (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::electric_rho is not implemented!");
}

template <class System> double CartesianField<System>::electric_phi (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::electric_phi is not implemented!");
}

template <class System> double CartesianField<System>::electric_r (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::electric_r is not implemented!");
}

template <class System> double CartesianField<System>::electric_theta (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::electric_theta is not implemented!");
}

template <class System> double CartesianField<System>::magnetic_rho (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::magnetic_rho is not implemented!");
}

template <class System> double CartesianField<System>::magnetic_phi (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::magnetic_phi is not implemented!");
}

template <class System> double CartesianField<System>::magnetic_r (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::magnetic_r is not implemented!");
}

template <class System> double CartesianField<System>::magnetic_theta (const Point::SpaceTime<System>&) const
{
    throw std::logic_error("CartesianField::magnetic_theta is not implemented!");
}

template <class System> double CartesianField<System>::energy   (const System&) const
{
    throw std::logic_error("CartesianField::energy is not implemented!");
}

template <class System> double CartesianField<System>::energy_e (const System&) const
{
    throw std::logic_error("CartesianField::energy_e is not implemented!");
}

template <class System> double CartesianField<System>::energy_h (const System&) const
{
    throw std::logic_error("CartesianField::energy_h is not implemented!");
}
