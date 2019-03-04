//
//  plot_energy_distribution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <string> // string() to_string()
#include <utility> // swap()
#include <iostream> // cout
using namespace std;

double R  = 1; // disk radius
double A0 = 1; // max current magnitude

double MU  = 1; // relative magnetic permatiity
double EPS = 1; // relative dielectric pirmativity

string MODULE = "module/uniform_disk/libuniform_disk.dylib"; // module path
int SUBMODULE = 0; // submodule index (meandr monocycle)
int THREAD_NUM = 5; // number of threads for calculation

auto str_of = [] (double val) { return to_string(val).substr(0,4); };

void plot_energy_distribution (double tau0, bool swap_axis = false)
{	
	ModuleManager mng = ModuleManager(NULL);
	mng.load_module(MODULE, R, A0, tau0, EPS, MU);
	// LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	// LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;

	// FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	// free_shape->set_time_depth([tau0] (double vt) {return Function::gauss(vt,tau0);});
	// LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);

	Manager<0>* thead_core = new Manager<0>(THREAD_NUM, NULL);
	thead_core->progress_bar(true);

	double x = 0;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= 5; z += 0.05) {
			double rho = sqrt(x*x + y*y);
			double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau0 + sqrt((rho+R)*(rho+R) + z*z) + 0.01;
			if (swap_axis) thead_core->add_argument( {y,x,z,from,to} );
			else thead_core->add_argument( {x,y,z,from,to} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = make_pair(property, linear);
	thead_core->call({function});
	auto data = thead_core->get_value();

	string script = "W_" + str_of(tau0) + (swap_axis?"_xz" :"_yz") + ".gnp";
	GnuPlot* plot = new GnuPlot(script);
	plot->set_colormap(Colormap::parula);
	plot->plot_colormap(data, !swap_axis, 2);

	delete linear, thead_core, plot;
}

int main ()
{
    plot_energy_distribution(2.0, false);
	plot_energy_distribution(2.0, true );
    return 0;
}
