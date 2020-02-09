//
//  space_point.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 21.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#define STRLEN 5

#include "phys_math.hpp"

#include <cmath>
#include <vector>
#include <string>
#include <initializer_list>

namespace Point {

    struct Polar;
    struct Spherical;
    struct Cylindrical;
    struct Cartesian3D;
    struct Cartesian2D;

    struct System : public std::vector<double> {
        using std::vector<double>::vector;
        virtual ~System () = default;
        virtual double radius () const { throw std::logic_error("Not implemented in Point::System"); }
        std::string to_str() const { std::string s; for (auto i : *this) s+=std::to_string(i).substr(0,STRLEN); return s; }
    private:
        using vector<double>::erase;
        using vector<double>::insert;
        using vector<double>::emplace;
        using vector<double>::pop_back;
        using vector<double>::push_back;
        using vector<double>::emplace_back;
    };

    template <class SystemImp> struct SpaceTime : public SystemImp {
        SpaceTime () = default;
        ~SpaceTime () = default;
        SpaceTime (const SystemImp& base) : SystemImp(base) { }
        SpaceTime (double ct, const SystemImp& base) : SystemImp(base), ctime(ct) { }
        SpaceTime (const std::initializer_list<double>::iterator& from, const std::initializer_list<double>::iterator& to) : SystemImp(from+1,to), ctime(*from)
        {
            if (std::distance(from,to) != 1+this->size()) throw std::logic_error("Number of arguments is not legal!");
            // TODO: bug!!! case never works
        }
        SpaceTime (const std::initializer_list<double> il) : SystemImp(il.begin()+1,il.end()), ctime(*il.begin()) 
        { if (il.size() != 1+this->size()) throw std::logic_error("Number of arguments is not legal!"); }
        double ct () const { return this->ctime; }
        double& ct () { return this->ctime; }
        double sqrt_vt2_z2 () const;
        bool casuality () { return this->ct() <= this->radius(); }
        std::string to_str() const { std::string s = std::to_string(ct()).substr(0,STRLEN); return s.append(SystemImp::to_str()); }
        template <class InputBase> static SpaceTime<SystemImp> convert (const SpaceTime<InputBase>& pt);
    private:
        double ctime{};
    };

    template <class SystemImp> struct ModalSpaceTime : public SpaceTime<SystemImp> {
        ModalSpaceTime () = default;
        ~ModalSpaceTime () = default;
        ModalSpaceTime (const SystemImp& base) : SpaceTime<SystemImp>(base) { }
        ModalSpaceTime (double m, double nu, const SpaceTime<SystemImp>& base) : SpaceTime<SystemImp>(base), _m(m), _nu(nu) { }
        ModalSpaceTime (const std::initializer_list<double> il) : SpaceTime<SystemImp>(il.begin()+2,il.end()), _m(*(il.begin())), _nu(*(il.begin()+1)) { }
        double m() const { return this->_m; }
        double& m() { return this->_m; }
        double nu() const { return this->_nu; }
        double& nu() { return this->_nu; }
    private:
        double _m;
        double _nu;
    };

    struct OneDim : public System {
        OneDim () : System(1,0) {}
        virtual ~OneDim () = default;
        OneDim (double val) : System(1,0) { at(0) = val; }
        double radius() const override { return at(0); }
        double& z() { return at(0); }
        double z() const { return at(0); }
    protected:
        using System::System;
    };

    struct Cartesian2D : public System {
        Cartesian2D () : System(2,0) {}
        Cartesian2D (double x, double y) : System(2,0) { at(0) = x; at(1) = y; }
        virtual ~Cartesian2D () = default;
        double radius () const override { return std::sqrt(x()*x() + y()*y()); }
        double& x() { return at(0); }
        double& y() { return at(1); }
        double x() const { return at(0); }
        double y() const { return at(1); }

        static Cartesian2D convert (const Polar& point);
        static Cartesian2D convert (const Spherical& point);
        static Cartesian2D convert (const Cartesian2D& point);
        static Cartesian2D convert (const Cartesian3D& point);
        static Cartesian2D convert (const Cylindrical& point);
    protected:
        using System::System;
    };

    struct Cartesian3D : public System {
        Cartesian3D () : System(3,0) {}
        Cartesian3D (double x, double y, double z) : System(3,0) { at(0) = x; at(1) = y; at(2) = z; }
        virtual ~Cartesian3D () = default;
        double radius () const override { return std::sqrt(x()*x() + y()*y() + z()*z()); }
        double& x() { return at(0); }
        double& y() { return at(1); }
        double& z() { return at(2); }
        double x() const { return at(0); }
        double y() const { return at(1); }
        double z() const { return at(2); }

        static Cartesian3D convert (const Polar& point);
        static Cartesian3D convert (const Spherical& point);
        static Cartesian3D convert (const Cartesian2D& point);
        static Cartesian3D convert (const Cartesian3D& point);
        static Cartesian3D convert (const Cylindrical& point);
    protected:
        using System::System;
    };

    struct Polar : public System {
        Polar () : System(2,0) {}
        Polar (double rho, double phi) : System(2,0) { at(0) = rho; at(1) = phi; }
        virtual ~Polar () = default;
        double radius () const override { return rho(); }
        double& phi() { return at(0); }
        double& rho() { return at(1); }
        double phi() const { return at(0); }
        double rho() const { return at(1); }

        static Polar convert (const Polar& point);
        static Polar convert (const Spherical& point);
        static Polar convert (const Cartesian2D& point);
        static Polar convert (const Cartesian3D& point);
        static Polar convert (const Cylindrical& point);
    protected:
        using System::System;
    };

    struct Cylindrical : public System {
        Cylindrical () : System(3,0) {}
        Cylindrical (double rho, double phi, double z) : System(3,0) { at(0) = rho; at(1) = phi; at(2) = z; }
        virtual ~Cylindrical () = default;
        double radius () const override { return std::sqrt(rho()*rho() + z()*z()); }
        double& rho() { return at(0); }
        double& phi() { return at(1); }
        double& z()   { return at(2); }
        double rho() const { return at(0); }
        double phi() const { return at(1); }
        double z()   const { return at(2); }

        static Cylindrical convert (const Polar& point);
        static Cylindrical convert (const Spherical& point);
        static Cylindrical convert (const Cartesian2D& point);
        static Cylindrical convert (const Cartesian3D& point);
        static Cylindrical convert (const Cylindrical& point);
    protected:
        using System::System;
    };

    struct Spherical : public System {
        Spherical() : System(3,0) {}
        Spherical (double r, double phi, double theta) : System(3,0) { at(0) = r; at(1) = phi; at(2) = theta; }
        virtual ~Spherical () = default;
        double radius () const override { return r(); }
        double& r()     { return at(0); }
        double& phi()   { return at(1); }
        double& theta() { return at(2); }
        double r()     const { return at(0); }
        double phi()   const { return at(1); }
        double theta() const { return at(2); }

        static Spherical convert (const Polar& point);
        static Spherical convert (const Spherical& point);
        static Spherical convert (const Cartesian2D& point);
        static Spherical convert (const Cartesian3D& point);
        static Spherical convert (const Cylindrical& point);
    protected:
        using System::System;
    };

    bool equals (const System& pt1, const System& pt2, double eps = DBL_MIN);
    template <class System> bool equals (const SpaceTime<System>& pt1, const SpaceTime<System>& pt2, double eps = DBL_MIN);
}

template <class SystemImp> double Point::SpaceTime<SystemImp>::sqrt_vt2_z2 () const
{
    return std::sqrt(ct()*ct()-this->z()*this->z());
}

template <class System> bool Point::equals (const Point::SpaceTime<System>& pt1, const Point::SpaceTime<System>& pt2, double eps)
{
    bool intime = Math::compare(pt1.ct(), pt2.ct(), eps);
    bool inspace = Point::equals(static_cast<System>(pt1),static_cast<System>(pt2),eps);
    return intime && inspace;
}

template <class System> template <class InputBase>
Point::SpaceTime<System> Point::SpaceTime<System>::convert (const Point::SpaceTime<InputBase>& pt)
{
    System base = System::convert(static_cast<InputBase>(pt));
    Point::SpaceTime<System> res{base};
    res.ct() = pt.ct();
    return res;
}
