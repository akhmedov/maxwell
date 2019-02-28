//
//  simpson2d_line.cpp
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
    static bool simpson2d_line ();
};

bool Test::simpson2d_line ()
{
	double a = 2;
	auto func = [] (double x, double y) { return x * y; };

	Simpson2D_line integral = Simpson2D_line();
	integral.first_limit(0, 400, a);
	auto min_y  = [] (double/*x*/) { return 0;     };
	auto term_y = [] (double x) { return (std::size_t) 200*x; };
	auto max_y  = [] (double x) { return x;     };
	integral.second_limit(min_y, term_y, max_y);

	double res = integral.value(func);
	double error = 100 * abs(res-2) / 2;
	return error < 1 ? true : false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::simpson2d_line()" << left;
    cout << (Test::simpson2d_line() ? "PASSED" : "FAILED") << endl;
}
