//
//  plot_energy_distribution.cpp
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

void plot_energy_distribution (double tau, double max_z)
{
	auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

	double R = 1, A0 = 1;
	double eps_r = 1, mu_r = 1;

	double x = 0;

	PlotTest::global_conf->field_component(FieldComponent::W);
	PlotTest::global_conf->impulse_shape(ImpulseShape::gauss);
	PlotTest::global_conf->duration(tau);

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::sinc(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

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
	auto function = std::make_pair(property, linear);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	GnuPlot* plot = new GnuPlot( str_of(tau) + "_" + str_of(max_z) + ".gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("y, m");
	plot->set_oy_label("z, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 1, 2);
}

int main ()
{
    global_conf = new Config();
	global_conf->path_gnuplot_binary("gnuplot/bin/gnuplot");
    plot_energy_distribution(0,2);
    return 0;
}
