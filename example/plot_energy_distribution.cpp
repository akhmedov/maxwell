//
//  plot_energy_distribution.cpp
//  example.interface.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

void plot_energy_distribution (double tau, double max_z)
{
	auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

	double R = 1, A0 = 1;
	double eps_r = 1, mu_r = 1;

	double x = 0;
	
	ModuleManager mng = ModuleManager();
	mng.load_module("module/uniform_disk/libuniform_disk.dylib");

	LinearCurrent* source = mng.get_module(mng.get_loaded()[0]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[0]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[0]).field;

	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);

	Manager<0>* thead_core = new Manager<0>(6, NULL);
	thead_core->progress_bar(true);

	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= max_z; z += 0.05) {
			double rho = std::sqrt(x*x + y*y);
			double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau + std::sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to+0.01} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = std::make_pair(property, duhamel);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	GnuPlot* plot = new GnuPlot( str_of(tau) + "_" + str_of(max_z) + ".gnp" );
	plot->set_gnuplot_bin("gnuplot/bin/gnuplot");
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("y, m");
	plot->set_oy_label("z, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 1, 2);
}

int main ()
{
    plot_energy_distribution(2,5);
    return 0;
}
