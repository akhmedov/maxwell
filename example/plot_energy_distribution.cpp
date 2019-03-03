//
//  plot_energy_distribution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

double R  = 1; // disk radius
double A0 = 1; // max current magnitude

double MU  = 1; // relative magnetic permatiity
double EPS = 1; // relative dielectric pirmativity

string MODULE = "module/uniform_disk/libuniform_disk.dylib"; // module path
int SUBMODULE = 0; // submodule index (meandr monocycle)
int THREAD_NUM = 5; // number of threads for calculation

auto str_of = [] (double val) { return to_string(val).substr(0,5); };

void plot_energy_distribution (double tau0, double max_z)
{	
	ModuleManager mng = ModuleManager(NULL);
	mng.load_module(MODULE, R, A0, tau0, EPS, MU);
	
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;

	// FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	// free_shape->set_time_depth([tau0] (double vt) {return Function::gauss(vt,tau0);});
	// LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);

	Manager<0>* thead_core = new Manager<0>(THREAD_NUM, NULL);
	thead_core->progress_bar(true);

	double y = 0;
	for (double x = -2*R; x <= 2*R; x += 0.05) {
		for (double z = 0; z <= 5; z += 0.05) {
			double rho = sqrt(x*x + y*y);
			double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
			double to = tau0 + sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = make_pair(property, linear);
	thead_core->call({function});
	vector<vector<double>> data = thead_core->get_value();

	GnuPlot* plot = new GnuPlot( "W_" + str_of(tau0) + ".gnp" );
	plot->set_gnuplot_bin("gnuplot/bin/gnuplot");
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("y, R");
	plot->set_oy_label("z, R");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 0, 2);
}

int main ()
{
    plot_energy_distribution(1,5);
    return 0;
}
