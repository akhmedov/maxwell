//
//  simpson2d.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"
using namespace std;

struct Test : Integral {
    static bool simpson2d ();
};

bool Test::simpson2d ()
{
	double radius = 1;

	auto func = [] (double rho, double/*phi*/) {
		return rho;
	};

	double etalon = M_PI * radius * radius;

	vector<tuple<double,size_t,double>> limits;
	limits.push_back( make_tuple(0,300,radius) );
	limits.push_back( make_tuple(0,300,M_PI_2) );

	Simpson2D integral = Simpson2D(limits);
	double square = 4 * integral.value(func);
	double error = 100 * abs(square-etalon) / etalon;

	return error < 1 ? true : false;
}

int main ()
{
	return !Test::simpson2d();
}
