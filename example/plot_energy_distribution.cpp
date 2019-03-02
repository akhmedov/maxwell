//
//  plot_energy_distribution.cpp
//  example.interface.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#define MODULE "module/uniform_disk/libuniform_disk.dylib"

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

void plot_energy_distribution (double tau0, double max_z)
{
	double R = 1, A0 = 1, eps_r = 1, mu_r = 1;
	
	ModuleManager mng = ModuleManager();
	mng.load_module(MODULE, NULL, R, A0, tau0, eps_r, mu_r);
	int submodule = 0;

	LinearCurrent* source = mng.get_module(mng.get_loaded()[submodule]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[submodule]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[submodule]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[submodule] << endl;

	// FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	// free_shape->set_time_depth([tau0] (double vt) {return Function::gauss(vt,tau0);});
	// LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);

	Manager<0>* thead_core = new Manager<0>(6, NULL);
	thead_core->progress_bar(true);

	double x = 0;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= max_z; z += 0.05) {
			double rho = std::sqrt(x*x + y*y);
			double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau0 + std::sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to+0.01} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = std::make_pair(property, linear);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	GnuPlot* plot = new GnuPlot( str_of(tau0) + "_" + str_of(max_z) + ".gnp" );
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
