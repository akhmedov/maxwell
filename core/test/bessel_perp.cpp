//
//  bessel_perp.cpp
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
    static bool bessel_perp ();
};

bool Test::bessel_perp ()
{
	auto f = [] (double x) { return j0(x); } ;
	auto f_perp = [] (double x) { return -j1(x); } ;

	size_t iterator = 0;
	double error = 0;
	for (double arg = -5; arg < 5; arg += 10e-5) {
		error = abs( f_perp(arg) - Math::derivat4(f,arg) );
		iterator++;
	}

	if (error/iterator > 10e-3) return false;
	else return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::bessel_perp()" << left;
    cout << (Test::bessel_perp() ? "PASSED" : "FAILED") << endl;
}
