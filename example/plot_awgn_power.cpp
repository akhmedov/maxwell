//
//  awgn_power.cpp
//  example.maxwell
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

void awgn_power (const vector<int>& samples)
{
	vector<vector<double>> data;
	vector<double> line;

	for  (double sigma = 0.1; sigma < 10; sigma += 0.1) {
		line = {sigma};
		WhiteGaussian* noise = new WhiteGaussian(0,sigma);
		for (auto i : samples) line.push_back(noise->power(i));
		data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( "awgn_power.gnp" );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("sigma");
	plot->set_oy_label("Pn, V*V");

	vector<string> title;
	for (auto i : samples) title.push_back("N = " + to_string(i));
	plot->plot_multi(data, title);
}

int main ()
{
    vector<int> samples{10,100,500};
    awgn_power(samples);
    return 0;
}
