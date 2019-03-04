//
//  binom_prod.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"

using namespace std;

struct Test : Math {
    static bool binom_prod ();
};

bool Test::binom_prod ()
{
	for (size_t n = 1; n <= 30; n++) {
		for (size_t m = 0; m <= n; m++) {
			size_t C1 = Math::binomial(n+1, m);
			size_t C2 = Math::binomial(n, n-m);
			size_t C12 = Math::binom_prod(n, m);
			if (C12 != C1*C2) return false;
		}
	}

	return true;
}

int main ()
{
	return !Test::binom_prod();
}