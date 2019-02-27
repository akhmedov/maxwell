//
//  monte_carlo_integral.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"

#include <iomanip>
#include <iostream>
using namespace std;

struct Test : Integral {
    static bool monte_carlo_integral ();
};

bool Test::monte_carlo_integral ()
{
	double radius = 1;

	auto func = [] (double rho, double phi, double theta) {
		UNUSED(phi);
		return rho * rho * sin(theta);
	};

	auto int_func = [radius] () {
		double R3 = radius * radius * radius;
		return 4 * M_PI * R3 / 3;
	};

	vector<pair<double,double>> limits;
	limits.push_back( make_pair(0,radius) );
	limits.push_back( make_pair(0,M_PI_2) );
	limits.push_back( make_pair(0,M_PI_2) );

	MonteCarlo integral = MonteCarlo (10e8, limits );
	double volume = 8 * integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 1 ? true : false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::monte_carlo_integral()" << left;
    // cout << (Test::monte_carlo_integral() ? "PASSED" : "FAILED") << endl;
	cout << "FAILED" << endl;
}