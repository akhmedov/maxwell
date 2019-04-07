//
//  space_point.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 21.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "space_point.hpp"

bool Point::equals (const Point::System& pt1, const Point::System& pt2, double eps)
{
    if (pt1.size() != pt2.size()) return false;
    for (auto i = 0u; i < pt1.size(); i++)
        if (!Math::compare(pt1[i], pt2[i], eps))
            return false;
    return true;
}

// Converters to Cartesian2D

Point::Cartesian2D Point::Cartesian2D::convert (const Point::Polar& pt)
{
    double x = pt.rho() * std::cos(pt.phi());
    double y = pt.rho() * std::sin(pt.phi());
    return Point::Cartesian2D(x, y);
}

Point::Cartesian2D Point::Cartesian2D::convert (const Point::Spherical& pt)
{
    double x = pt.r() * std::sin(pt.theta()) * std::cos(pt.phi());
    double y = pt.r() * std::sin(pt.theta()) * std::sin(pt.phi());
    return Point::Cartesian2D(x,y);
}

Point::Cartesian2D Point::Cartesian2D::convert (const Point::Cartesian2D& pt)
{
    return pt;
}

Point::Cartesian2D Point::Cartesian2D::convert (const Point::Cartesian3D& pt)
{
    return Point::Cartesian2D(pt.x(), pt.y());
}

Point::Cartesian2D Point::Cartesian2D::convert (const Point::Cylindrical& pt)
{
    double x = pt.rho() * std::cos(pt.phi());
    double y = pt.rho() * std::sin(pt.phi());
    return Point::Cartesian2D(x, y);
}

// Converters to Cartesian3D

Point::Cartesian3D Point::Cartesian3D::convert (const Point::Cartesian3D& pt)
{
    return pt;
}

Point::Cartesian3D Point::Cartesian3D::convert (const Point::Cylindrical& pt)
{
    double x = pt.rho() * std::cos(pt.phi());
    double y = pt.rho() * std::sin(pt.phi());
    double z = pt.z();
    return Point::Cartesian3D(x,y,z);
}

Point::Cartesian3D Point::Cartesian3D::convert (const Point::Spherical& pt)
{
    double x = pt.r() * std::sin(pt.theta()) * std::cos(pt.phi());
    double y = pt.r() * std::sin(pt.theta()) * std::sin(pt.phi());
    double z = pt.r() * std::cos(pt.theta());
    return Point::Cartesian3D(x,y,z);
}

Point::Cartesian3D Point::Cartesian3D::convert (const Point::Polar& pt)
{
    double x = pt.rho() * std::cos(pt.phi());
    double y = pt.rho() * std::sin(pt.phi());
    return Point::Cartesian3D(x,y,0);
}

Point::Cartesian3D Point::Cartesian3D::convert (const Point::Cartesian2D& pt)
{
    return Point::Cartesian3D(pt.x(), pt.y(), 0);
}

// Converters to Cylindrical

Point::Cylindrical Point::Cylindrical::convert (const Point::Cartesian3D& pt)
{
    double rho = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y());
    double phi = std::atan2(pt.y(),pt.x());
    double z = pt.z();
    return Point::Cylindrical(rho,phi,z);
}

Point::Cylindrical Point::Cylindrical::convert (const Point::Spherical& pt)
{
    double rho = pt.r() * std::sin(pt.theta());
    double phi = pt.phi();
    double z = pt.r() * std::cos(pt.theta());
    return Point::Cylindrical(rho,phi,z);
}

Point::Cylindrical Point::Cylindrical::convert (const Point::Cylindrical& pt)
{
    return pt;
}

Point::Cylindrical Point::Cylindrical::convert (const Point::Polar& pt)
{
    return Point::Cylindrical(pt.rho(),pt.phi(),0);
}

Point::Cylindrical Point::Cylindrical::convert (const Point::Cartesian2D& pt)
{
    double rho = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y());
    double phi = std::atan2(pt.y(),pt.x());
    return Point::Cylindrical(rho,phi,0);
}

// Converters to Spherical

Point::Spherical Point::Spherical::convert (const Point::Cartesian3D& pt)
{
    double r = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
    double phi = std::atan2(pt.y(),pt.x());
    double theta = std::acos( pt.z() / r );
    return Point::Spherical(r,phi,theta);
}

Point::Spherical Point::Spherical::convert (const Point::Cylindrical& pt)
{
    double r = std::sqrt(pt.rho() * pt.rho() + pt.z() * pt.z());
    double phi = pt.phi();
    double theta = std::atan2( pt.rho(), pt.z() );
    return Point::Spherical(r,phi,theta);
}

Point::Spherical Point::Spherical::convert (const Point::Spherical& pt)
{
    return pt;
}

Point::Spherical Point::Spherical::convert (const Point::Cartesian2D& pt)
{
    double r = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y());
    double phi = std::atan2(pt.y(),pt.x());
    double theta = std::acos( 0 );
    return Point::Spherical(r,phi,theta);
}

Point::Spherical Point::Spherical::convert (const Point::Polar& pt)
{
    double r = std::sqrt(pt.rho() * pt.rho());
    double phi = pt.phi();
    double theta = std::atan2( pt.rho(), 0 );
    return Point::Spherical(r,phi,theta);
}