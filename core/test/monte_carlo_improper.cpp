//
//  monte_carlo_improper.cpp
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
    static bool monte_carlo_improper ();
};

bool Test::monte_carlo_improper ()
{
	double a = 1, b = 1;

	auto func = [a,b] (double x, double y, double z) {
		return exp(-y-z) * j0(a*sqrt(x*y)) * j0(b*sqrt(x*z));
	};

	auto int_func = [a,b] () {
		return a*a/4 + b*b/4;
	};

	vector<pair<double,double>> limits;
	limits.push_back(make_pair(0,10e3));
	limits.push_back(make_pair(0,10e3));
	limits.push_back(make_pair(0,10e3));

	MonteCarlo integral = MonteCarlo(10e8, limits);
	double volume = integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 10 ? true : false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::monte_carlo_improper()" << left;
    // cout << (Test::monte_carlo_improper() ? "PASSED" : "FAILED") << endl;
	cout << "FAILED" << endl;
}
