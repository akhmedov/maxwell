//
//  awgn_power.cpp
//  example.interface.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "noise.hpp"
#include "config.hpp"
#include "gnu_plot.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

Config* global_conf;

void awgn_power (const vector<int>& samples)
{
	vector<vector<double>> data;
	vector<double> line;

	for  (double sigma = 0.1; sigma < 10; sigma += 0.1) {
		line = {sigma};
		AdditiveWhiteGaussian* noise = new AdditiveWhiteGaussian(0,sigma);
		for (auto i : samples) {
			double Pn = noise->power(0,0,0,i);
			line.push_back(Pn);
		}
		data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( "power.gnp" );
	plot->set_gnuplot_bin( global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("sigma");
	plot->set_oy_label("Pn, V*V");
	plot->grid_on(false);
	plot->cage_on();
	vector<string> title;
	for (auto i : samples) title.push_back("N = " + to_string(i));
	plot->plot_multi(data, title);
}

int main ()
{
    global_conf = new Config();
	global_conf->path_gnuplot_binary("gnuplot/bin/gnuplot");
    vector<int> samples = {10,100,500};
    awgn_power(samples);
    return 0;
}
