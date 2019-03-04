//
//  simpson3d.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"
using namespace std;

struct Test : Integral {
    static bool simpson3d ();
};

bool Test::simpson3d ()
{
	double radius = 1;

	auto func = [] (double rho, double/*phi*/, double theta) {
		return rho * rho * sin(theta);
	};

	auto int_func = [radius] () {
		double R3 = radius * radius * radius;
		return 4 * M_PI * R3 / 3;
	};

	vector<tuple<double,size_t,double>> limits;
	limits.push_back( make_tuple(0,200,radius) );
	limits.push_back( make_tuple(0,200,M_PI_2) );
	limits.push_back( make_tuple(0,200,M_PI_2) );

	Simpson3D integral = Simpson3D(limits);
	double volume = 8 * integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 1 ? true : false;
}

int main ()
{
	return !Test::simpson3d();
}
