//
//  hardcode_dataset.cpp
//  example.interface.maxwell
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

bool hardcode_dataset ()
{
	double tau = 0.5;
	double duty_cycle = 0.5;
	double vt_step = 0.01;

	std::pair<double,double> rho = std::make_pair(0,7);
	std::pair<double,double> phi = std::make_pair(0,90);
	std::pair<double,double> z = std::make_pair(0,20);

	std::vector<std::function<double(double)>> domain = {
		[tau] (double vt) { return   Function::sinc  (vt,tau); },
		[tau] (double vt) { return   Function::gauss (vt,tau); }
	};

	serial::dataset ds = serial::dataset(domain, tau, NULL, duty_cycle, vt_step);
	
	ds.set_char(0, 0, 2, 1);
	ds.set_char(0, 0, 2, 2);
	ds.set_char(0, 0, 2, 2);
	ds.set_char(0, 0, 2, 1);

	json js = serial::json_from(ds,rho,phi,z);
	serial::serialize("dataset.json", js);
	return true;
}

int main ()
{
    hardcode_dataset();
    return 0;
}
