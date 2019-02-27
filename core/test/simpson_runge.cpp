//
//  simpson_runge.cpp
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
    static bool simpson_runge ();
};

bool Test::simpson_runge ()
{
	double accuracy = 1; // %

	double from = 0;
	double to = M_PI_2;

	auto f = [] (double x) {
		return std::sin(x);
	};

	double I = SimpsonRunge(1, 1).value(from, to, f);
	return 100*std::abs(I-1) < accuracy ? true : false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::simpson_runge()" << left;
    cout << (Test::simpson_runge() ? "PASSED" : "FAILED") << endl;
}
