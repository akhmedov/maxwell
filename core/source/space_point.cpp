//
//  space_point.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 21.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "space_point.hpp"

Point::Cartesian2D::Cartesian2D(double x, double y)
{ 
    this->at(0) = x;
    this->at(1) = y;
}

double Point::Cartesian2D::radius () const
{
    return std::sqrt(x()*x() + y()*y());
}

Point::Cartesian3D::Cartesian3D(double x, double y, double z)
{
    this->at(0) = x;
    this->at(1) = y;
    this->at(2) = z;
}

double Point::Cartesian3D::radius () const
{
    return std::sqrt(x()*x() + y()*y() + z()*z());
}

Point::Polar::Polar(double rho, double phi)
{
    this->at(0) = rho;
    this->at(1) = phi;
}

double Point::Polar::radius () const
{
    return std::abs(rho());
}

Point::Cylindrical::Cylindrical(double rho, double phi, double z)
{
    this->at(0) = rho;
    this->at(1) = phi;
    this->at(2) = z;
}

double Point::Cylindrical::radius () const
{
    return std::sqrt(rho() * rho() + z() * z());
}

Point::Spherical::Spherical(double r, double phi, double theta)
{
    this->at(0) = r;
    this->at(1) = phi;
    this->at(2) = theta;
}

double Point::Spherical::radius () const
{
    return std::abs(r());
}

Point::Cartesian3D Point::cartesian3d (const Point::Cylindrical& pt)
{
    Point::Cartesian3D res;
    res.x() = pt.rho() * std::cos(pt.phi());
    res.y() = pt.rho() * std::sin(pt.phi());
    res.z() = pt.z();
    return res;
}

Point::Cartesian3D Point::cartesian3d (const Point::Spherical& pt)
{
    Point::Cartesian3D res;
    res.x() = pt.r() * std::sin(pt.theta()) * std::cos(pt.phi());
    res.y() = pt.r() * std::sin(pt.theta()) * std::sin(pt.phi());
    res.z() = pt.r() * std::cos(pt.theta());
    return res;
}

Point::Cylindrical Point::cylindrical (const Point::Cartesian3D& pt)
{
    Point::Cylindrical res;
    res.rho() = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y());
    res.phi() = std::atan2(pt.y(),pt.x());
    res.z() = pt.z();
    return res;
}

Point::Cylindrical Point::cylindrical (const Point::Spherical& pt)
{
    Point::Cylindrical res;
    res.rho() = pt.r() * std::sin(pt.theta());
    res.phi() = pt.phi();
    res.z() = pt.r() * std::cos(pt.theta());
    return res;
}

Point::Spherical Point::spherical (const Point::Cartesian3D& pt)
{
    Point::Spherical res;
    res.r() = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
    res.phi() = std::atan2(pt.y(),pt.x());
    res.theta() = std::acos( pt.z() / res.r() );
    return res;
}

Point::Spherical Point::spherical (const Point::Cylindrical& pt)
{
    Point::Spherical res;
    res.r() = std::sqrt(pt.rho() * pt.rho() + pt.z() * pt.z());
    res.phi() = pt.phi();
    res.theta() = std::atan2( pt.rho(), pt.z() );
    return res;
}
