//
//  plot_energy_slyse.cpp
//  example.interface.maxwell
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

void plot_energy_slyse (double tau, double z)
{
	auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

	double R = 1, A0 = 1;
	double eps_r = 1, mu_r = 1;
	double range = z/2;

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::sinc(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

	Manager<0>* thead_core = new Manager<0>(6, NULL);
	thead_core->progress_bar(true);

	for (double x = -range; x <= range; x += 0.05) {
		for (double y = -range; y <= range; y += 0.05) {
			double rho = std::sqrt(x*x + y*y);
			double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau + std::sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to+0.01} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = std::make_pair(property, linear);
	// auto function = std::make_pair(property, duhamel);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	std::vector<double> max0 = data[0];
	for (auto point : data)
		if (point[3] > max0[3])
			max0 = point;

	// for (auto&& i : data) {
	// 	i[3] *= z*z / max0[3]; // norm W
	// 	i.erase(i.begin()+2,i.begin()+5); // erase z, from, to
	// }

	std::cout << "Wmax (" << str_of(max0[0]) << ',' << str_of(max0[1]) << ',' << str_of(max0[2]) << ") = " << str_of(max0[3]) << std::endl;

	/* std::vector<std::vector<double>> agumented;
	for (auto&& i : data) {
		double x = i[0], y = i[1], W = i[2];
		if (x > 1e-8) agumented.push_back( {-x,y,W} );
		if (y > 1e-8) agumented.push_back( {x,-y,W} );
		if (x > 1e-8 && y > 1e-8) agumented.push_back( {-x,-y,W} );
	}
	data.insert(std::end(data), std::begin(agumented), std::end(agumentat)); */

	GnuPlot* plot = new GnuPlot( str_of(tau) + "_" + str_of(z) + ".gnp" );
	plot->set_gnuplot_bin("gnuplot/bin/gnuplot");
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("x, m");
	plot->set_oy_label("y, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 0, 1);
}

int main ()
{
    plot_energy_slyse(0,2);
    return 0;
}
