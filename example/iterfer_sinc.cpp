//
//  iterfer_sinc.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"
#include "uniform_disk_current.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

void iterfer_sinc ()
{
	double R = 1, A0 = 1, tau = 0.5;
	double eps_r = 1, mu_r = 1;
	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);
	auto property = &AbstractField::energy_cart;
	auto compute = make_pair(property, (AbstractField*)duhamel);

	Manager<2>* thead_core = new Manager<2>(4, NULL);
	thead_core->progress_bar(true);
	for (double z = 8; z < 16; z += 0.1) {
		double rho = 0;
		double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
		double to = tau + sqrt((rho+R)*(rho+R) + z*z);
		thead_core->add_argument( {0,0,z,from,to} );
	}

	thead_core->call({compute});
	vector<vector<double>> data = thead_core->get_value();

	for (auto&& i : data) {
		i.erase(i.begin(), i.begin() + 2); // erase x y
		// i[1] *= i[0] * i[0];
	}

	GnuPlot* plot = new GnuPlot( "iterfer_sinc.gnp" );
	plot->set_gnuplot_bin("gnuplot/bin/gnuplot");
	plot->set_ox_label("z, m");
	plot->set_oy_label("W, V*V/m2");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(data);
	plot->call_gnuplot();
}

int main ()
{
    iterfer_sinc();
    return 0;
}
