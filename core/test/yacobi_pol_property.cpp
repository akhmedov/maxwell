//
//  yacobi_pol_property.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

struct Test : Math {
    static bool yacobi_pol_property ();
};

bool Test::yacobi_pol_property ()
{
	double eps = 10e-6;

	double table_value = Math::binomial(5,4);
	double value_basic = Math::yacobi_polinom_init(4, 1);
	double value_fast = Math::yacobi_polinom(4, 1).get_d();

	if (abs(table_value - value_basic) > eps) return false;
	if (abs(table_value - value_fast) > eps) return false;

	table_value = pow(-1,4) * Math::binomial(4, 4);
	value_basic = Math::yacobi_polinom_init(4, -1);
	value_fast = Math::yacobi_polinom(4, -1).get_d();

	if (abs(table_value - value_basic) > eps) return false;
	if (abs(table_value - value_fast) > eps) return false;

	size_t order = 4;
	size_t arg = 5;
	double res_basic = Math::yacobi_polinom_init(order, arg);
	double res_fast = Math::yacobi_polinom(order, arg).get_d();

	if (std::abs(res_basic - res_fast) > eps) return false;

	for (double arg = -1; arg < 1; arg += 0.01) {
		double res_fast = Math::yacobi_polinom(0, arg).get_d();
		if (std::abs(res_fast - 1) > eps) return false;
	}

	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::yacobi_pol_property()" << left;
    cout << (Test::yacobi_pol_property() ? "PASSED" : "FAILED") << endl;
    return 0;
}
