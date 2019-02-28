//
//  plot_energy_max.cpp
//  example.interface.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "uniform_disk_current.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;


void PlotTest::plot_energy_max ()
{	
	double R = 1, A0 = 1;

	auto field = [R,A0] (double tau) {
		double eps_r = 1, mu_r = 1;
		Homogeneous* medium = new Homogeneous(mu_r, eps_r);
		UniformPlainDisk* source = new UniformPlainDisk(R, A0);
		MissileField* linear = new MissileField(source, medium);
		FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
		free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
		LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);
		auto property = &AbstractField::energy_cart;
		return std::make_pair(property, duhamel);
	};

	Manager<2>* thead_core = new Manager<2>(4, NULL);
	thead_core->progress_bar(true);

	for (double z = 0.5; z <= 15; z += 0.1) {
		double rho = 0;
		double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
		double to = 1 + std::sqrt((rho+R)*(rho+R) + z*z);
		thead_core->add_argument( {0,0,z,from,to} );
	}

	thead_core->call({field(0.5), field(1)});
	std::vector<std::vector<double>> data = thead_core->get_value();
	
	double max0 = 0;
	for (auto i : data) 
		if (std::abs(i[2] - 15) < 1e-5) 
			max0 = i[3] * i[2] * i[2];

	for (auto&& i : data) {
		i[3] *= i[2] * i[2] / max0; // norm W
		i[4] *= i[2] * i[2] / max0; // norm W
		i.erase(i.begin(), i.begin() + 2); // erase x and y
	}

	GnuPlot* plot = new GnuPlot("plot_energy_max.gnp");
	plot->set_gnuplot_bin("gnuplot/bin/gnuplot");
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("z, m");
	plot->set_oy_label("W, V*V");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"tau = 0.5 vt", 
									  "tau = 1.0 vt"};
	plot->plot_multi(data, title);
}

int main ()
{
    plot_energy_max(0,2);
    return 0;
}
