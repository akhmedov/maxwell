//
//  monte_carlo_vector.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 26.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"
using namespace std;

struct Test : Integral {
    static bool monte_carlo_vector ();
};

bool Test::monte_carlo_vector ()
{
	vector<pair<double,double>> limits;
	limits.push_back(make_pair(0,1));
	limits.push_back(make_pair(0,1));

	MonteCarlo distr = MonteCarlo (3000, limits );
	std::valarray<double> val = distr.random_array();

	return false;
}

int main ()
{
	// return !Test::monte_carlo_vector();
	return 0;
}