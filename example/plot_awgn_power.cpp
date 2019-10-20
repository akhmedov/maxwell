//
//  awgn_power.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

#include <vector>
#include <iostream>
using namespace std;

void awgn_power (const vector<int>& samples)
{
	vector<string> title;
	vector<vector<double>> arg, fnc;

	for (auto i : samples) {
		title.push_back("N = " + to_string(i));
		vector<double> X, Y;
		for (double sigma = 0.1; sigma < 10; sigma += 0.1) {
			WhiteGaussian* noise = new WhiteGaussian(0,sigma);
			X.push_back(sigma);
			Y.push_back(noise->power(i));
		}
		arg.push_back(X);
		fnc.push_back(Y);
	}

	PyPlotManager plot = PyPlotManager("awgn_power.py");
	plot.set_ox_label("sigma");
	plot.set_oy_label("Pn, V*V");
	plot.set_colormap(ScriptManager::Colormap::grey);
	plot.plot2d(arg, fnc, title);
}

int main ()
{
    vector<int> samples{10,100,500};
    awgn_power(samples);
    return 0;
}
