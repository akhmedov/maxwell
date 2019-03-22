//
//  point_converter.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 21.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include <iostream>
#include "phys_math.hpp"
#include "space_point.hpp"

using namespace std;

static const double eps = 10e-9;

static const Point::SpaceTime<Point::Cartesian3D> cart{{10, M_LN10, M_PI, M_E}};
static const Point::SpaceTime<Point::Cylindrical> cyln{{10, M_LN10, M_PI, M_E}};
static const Point::SpaceTime<Point::Spherical>   sphe{{10, M_LN10, M_PI, M_E}};

template <class System> void print_point (const System& pt)
{
    for (auto i : pt) cout << i << ' ';
    cout << endl;
}

template <class System> void print_point (const Point::SpaceTime<System>& pt)
{
    cout << pt.ct() << ' ';
    print_point((System) pt);
}

int main ()
{
    Point::SpaceTime<Point::Space<3>> conv;

    conv = Point::cartesian3d(Point::cylindrical(cart));
    if (!Point::equals(conv,cart,eps)) return 1;

    conv = Point::cartesian3d(Point::spherical(cart));
    if (!Point::equals(conv,cart,eps)) return 1;

    conv = Point::cylindrical(Point::spherical(cyln));
    if (!Point::equals(conv,cyln,eps)) return 1;

    conv = Point::cylindrical(Point::cartesian3d(cyln));
    if (!Point::equals(conv,cyln,eps)) return 1;

    conv = Point::spherical(Point::cylindrical(sphe));
    if (!Point::equals(conv,sphe,eps)) return 1;

    conv = Point::spherical(Point::cartesian3d(sphe));
    if (!Point::equals(conv,sphe,eps)) return 1;
    
    return 0;
}
