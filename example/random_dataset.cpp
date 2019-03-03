//
//  random_dataset.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "dataset.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

bool random_dataset ()
{
	double radix = 3;
	double sigma = 5.0;
	double pulses = 60;
	std::string file_name = "dataset.json";
	std::pair<double,double> rho = std::make_pair(0,4);
	std::pair<double,double> phi = std::make_pair(0,20);
	std::pair<double,double> z = std::make_pair(4,6);
	serial::randomized_sequental (pulses, radix, sigma, file_name, rho, phi, z);
	return true;
}

int main ()
{
    random_dataset();
    return 0;
}
