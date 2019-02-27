//
//  binom_prod.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"

#include <iomanip>
#include <iostream>
using namespace std;

struct Test : Math {
    static bool yacobi_pol_compare ();
};

bool Test::yacobi_pol_compare ()
{
	double eps = 10e-4;

	for (size_t order = 1; order < 18; order++) {
		for (double arg = -3; arg < 3; arg += 0.01) {
			double res_basic = Math::yacobi_polinom_init(order, arg);
			double res_fast = Math::yacobi_polinom(order, arg).get_d();

			if (std::abs(res_basic - res_fast) > eps) return false;
		}		
	}

	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::yacobi_pol_compare()" << left;
    cout << (Test::yacobi_pol_compare() ? "PASSED" : "FAILED") << endl;
}
