//
//  invers_sqrt.cpp
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
    static bool invers_sqrt ();
};

bool Test::invers_sqrt ()
{
	double val = 2;
	double sqr = Math::inv_sqrt (val);
	if ( abs(sqr - 0.707) > 10e-5 ) return false;
	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::invers_sqrt()" << left;
    cout << (Test::invers_sqrt() ? "PASSED" : "FAILED") << endl;
}