//
//  hardcode_dataset.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "dataset.hpp"
#include "uniform_disk_current.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

AbstractField* updisk_field (const std::function<double(double)>& f) 
{
	double eps_r = 1, mu_r = 1;
	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(RADIUS, AMPLITUDE);
	MissileField* on = new MissileField(source, medium);
	auto property = &AbstractField::electric_x;
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source); 
	free_shape->set_time_depth(f);	
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, on, NULL);
	return (AbstractField*)duhamel;
}

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

	std::vector<AbstractField*> filed;
	for (auto i : domain) field.push_back(updisk_field(i));

	serial::dataset ds = serial::dataset(filed, tau, NULL, duty_cycle, vt_step);
	
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
