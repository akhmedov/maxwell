//
//  iterfer_sinc.cpp
//  example.interface.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "manager.hpp"
#include "function.hpp"

#include "config.hpp"
#include "gnu_plot.hpp"
#include "linear_duhamel.hpp"

#include "uniform_disk_current.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

Config* global_conf;

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

	GnuPlot* plot = new GnuPlot( global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( global_conf->path_gnuplot_binary() );
	plot->set_ox_label("z, m");
	plot->set_oy_label("W, V*V/m2");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(data);
	plot->call_gnuplot();
}

int main ()
{
    global_conf = new Config();
	global_conf->path_gnuplot_binary("gnuplot/bin/gnuplot");
    iterfer_sinc();
    return 0;
}
