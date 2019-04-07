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
using namespace Point;

static const double eps = 10e-9;

static const SpaceTime<Cartesian3D> CART{10, M_LN10, M_PI, M_E};
static const SpaceTime<Cylindrical> CYLN{10, M_LN10, M_PI, M_E};
static const SpaceTime<Spherical>   SPHE{10, M_LN10, M_PI, M_E};

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
    SpaceTime<Cartesian3D> cart;
    cart = SpaceTime<Cartesian3D>::convert(SpaceTime<Cylindrical>::convert(CART));
    if (!equals(CART,cart,eps)) return 1;
    cart = SpaceTime<Cartesian3D>::convert(SpaceTime<Spherical>::convert(CART));
    if (!equals(CART,cart,eps)) return 1;

    SpaceTime<Cylindrical> cyln;
    cyln = SpaceTime<Cylindrical>::convert(SpaceTime<Spherical>::convert(CYLN));
    if (!equals(CYLN,cyln,eps)) return 1;
    cyln = SpaceTime<Cylindrical>::convert(SpaceTime<Cartesian3D>::convert(CYLN));
    if (!equals(CYLN,cyln,eps)) return 1;

    SpaceTime<Spherical>   sphe;
    sphe = SpaceTime<Spherical>::convert(SpaceTime<Cylindrical>::convert(SPHE));
    if (!equals(SPHE,sphe,eps)) return 1;
    sphe = SpaceTime<Spherical>::convert(SpaceTime<Cartesian3D>::convert(SPHE));
    if (!equals(SPHE,sphe,eps)) return 1;
    
    return 0;
}
