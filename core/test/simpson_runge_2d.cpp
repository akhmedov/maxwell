//
//  simpson_runge_2d.cpp
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
    static bool simpson_runge_2d ();
};

bool Test::simpson_runge_2d ()
{
	size_t init_units = 1;
	double error = 2; /* % */

	auto f = [] (double x, double y) {
		return 4 * x * y;
	};

	SimpsonRunge integral = SimpsonRunge(init_units, error);

	double I = integral.value(0, 1,
		[f, &integral] (double x) {
			return integral.value(0, 1, [f,x] (double y) {
				return f(x,y); 
			} );
		} );

	return 100 * abs(I - 1) < 2*error ? true : false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::simpson_runge_2d()" << left;
    cout << (Test::simpson_runge_2d() ? "PASSED" : "FAILED") << endl;
}
