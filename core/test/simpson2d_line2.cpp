//
//  simpson2d_line2.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"
using namespace std;

struct Test : Integral {
    static bool simpson2d_line2 ();
};

bool Test::simpson2d_line2 ()
{
	double radius = 1;

	auto func = [] (double rho, double/*phi*/) {
		return rho;
	};

	Simpson2D_line integral = Simpson2D_line();
	integral.first_limit(0, 300, radius);
	auto min_y  = [] (double/*x*/) { return 0;      };
	auto term_y = [] (double/*x*/) { return 300;    };
	auto max_y  = [] (double/*x*/) { return M_PI_2; };
	integral.second_limit(min_y, term_y, max_y);

	double res = 4 * integral.value(func);
	double error = 100 * abs(res-M_PI) / M_PI;
	return error < 1 ? true : false;
}

int main ()
{
	return !Test::simpson2d_line2();
}
