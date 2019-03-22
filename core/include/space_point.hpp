//
//  space_point.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 21.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "phys_math.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace Point {

    template <std::size_t dim> struct Space : public std::array<double,dim> {
        using std::array<double,dim>::array;
        Space (const std::vector<double>& data);
        virtual double radius () const;
    };

    // TODO: implement converters
    struct SingleDim : public Space<1> {
        using Space<1>::Space;
        SingleDim(double z);
        double radius () const;
        double& z() { return at(0); }
    };

    template <class System> struct SpaceTime : public System {
        using System::System;
        // SpaceTime () = default;
        SpaceTime (const System& base);
        SpaceTime (const std::vector<double>& data);
        bool casuality ();
        double& ct() { return time; }
        double ct() const { return time; }
    private:
        double time{};
    };

    // TODO: implement converters
    struct Cartesian2D : public Space<2> {
        using Space<2>::Space;
        Cartesian2D(double x, double y);
        double radius () const;
        double& x() { return at(0); }
        double& y() { return at(1); }
        double x() const { return at(0); }
        double y() const { return at(1); }
    };

    struct Cartesian3D : public Space<3> {
        using Space<3>::Space;
        Cartesian3D(double x, double y, double z);
        double radius () const;
        double& x() { return at(0); }
        double& y() { return at(1); }
        double& z() { return at(2); }
        double x() const { return at(0); }
        double y() const { return at(1); }
        double z() const { return at(2); }
    };

    // TODO: implement converters
    struct Polar : public Space<2> {
        using Space<2>::Space;
        Polar(double rho, double phi);
        double radius () const;
        double& phi() { return at(0); }
        double& rho() { return at(1); }
        double phi() const { return at(0); }
        double rho() const { return at(1); }
    };

    struct Cylindrical : public Space<3> {
        using Space<3>::Space;
        Cylindrical(double rho, double phi, double z);
        double radius () const;
        double& rho() { return at(0); }
        double& phi() { return at(1); }
        double& z()   { return at(2); }
        double rho() const { return at(0); }
        double phi() const { return at(1); }
        double z()   const { return at(2); }
    };

    struct Spherical : public Space<3> {
        using Space<3>::Space;
        Spherical(double r, double phi, double theta);
        double radius () const;
        double& r()     { return at(0); }
        double& phi()   { return at(1); }
        double& theta() { return at(2); }
        double r()     const { return at(0); }
        double phi()   const { return at(1); }
        double theta() const { return at(2); }
    };

    Cartesian3D cartesian3d (const Cylindrical& pt);
    Cartesian3D cartesian3d (const Spherical& pt);
    template <class System> SpaceTime<Cartesian3D> cartesian (const SpaceTime<System>& pt);

    Cylindrical cylindrical (const Cartesian3D& pt);
    Cylindrical cylindrical (const Spherical& pt);
    template <class System> SpaceTime<Cylindrical> cylindrical (const SpaceTime<System>& pt);

    Spherical spherical (const Cartesian3D& pt);
    Spherical spherical (const Cylindrical& pt);
    template <class System> SpaceTime<Spherical> spherical (const SpaceTime<System>& pt);

    template <class System> bool equals (const SpaceTime<System>& pt1, const SpaceTime<System>& pt2, double eps = DBL_MIN);
    template <std::size_t dim> bool equals (const Space<dim>& pt1, const Space<dim>& pt2, double eps = DBL_MIN);
};

template <std::size_t dim> Point::Space<dim>::Space (const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i++)
        this->at(i) = data.at(i);
}

template <class System> Point::SpaceTime<System>::SpaceTime (const System& base)
{
    for (std::size_t i = 0; i < base.size(); i++)
        this->at(i) = base.at(i);
}

template <class System> Point::SpaceTime<System>::SpaceTime (const std::vector<double>& data)
{
    this->ct() = data[0];
    for (std::size_t i = 1; i < data.size(); i++) this->at(i-1) = data.at(i);
}

template <class System> bool Point::SpaceTime<System>::casuality ()
{
    return !(this->ct() < this->radius());
}

template <class System> bool Point::equals (const Point::SpaceTime<System>& pt1, const Point::SpaceTime<System>& pt2, double eps)
{
    bool intime = Math::compare(pt1.ct(), pt2.ct(), eps);
    bool inspace = Point::equals((System)pt1,(System)pt2,eps);
    return intime && inspace;
}

template <std::size_t dim> bool Point::equals (const Point::Space<dim>& pt1, const Point::Space<dim>& pt2, double eps)
{
    for (auto i = 0u; i < dim; i++)
        if (!Math::compare(pt1[i], pt2[i], eps))
            return false;
    return true;
}

template <class System> Point::SpaceTime<Point::Cartesian3D> Point::cartesian (const Point::SpaceTime<System>& pt)
{
    Point::Cartesian3D base = cartesian((System) pt);
    Point::SpaceTime<Point::Cartesian3D> res{base};
    res.ct() = pt.ct();
    return res;
}

template <class System> Point::SpaceTime<Point::Cylindrical> Point::cylindrical (const Point::SpaceTime<System>& pt)
{
    Point::Cylindrical base = cylindrical((System) pt);
    Point::SpaceTime<Point::Cylindrical> res{base};
    res.ct() = pt.ct();
    return res;
}

template <class System> Point::SpaceTime<Point::Spherical> Point::spherical (const Point::SpaceTime<System>& pt)
{
    Point::Spherical base = spherical((System) pt);
    Point::SpaceTime<Point::Spherical> res{base};
    res.ct() = pt.ct();
    return res;
}
