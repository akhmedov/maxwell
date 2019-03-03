//
//  snr_dataset.cpp
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

bool snr_dataset ()
{
	double radix = 3;
	double snr = 40;
	double pulses = 1e3;
	std::string file_name = "train.json";
	double rho = 0;
	double phi = 0;
	double z = 10;
	serial::same_snr(pulses, radix, snr, file_name, rho, phi, z);
	return true;
}

int main ()
{
    snr_dataset();
    return 0;
}
