//
//  plot_energy_distribution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"
#include "data.hpp"
#include "cluster.hpp"

#include <string> // string() to_string()
#include <utility> // swap()
#include <iostream> // cout
using namespace std;

double R  = 1; // disk radius
double A0 = 1; // max current magnitude

double MU  = 1; // relative magnetic permatiity
double EPS = 1; // relative dielectric pirmativity

auto str_of = [] (double val) { return to_string(val).substr(0,4); };

AbstractField* create_model (double tau0) 
{
	string MODULE = "module/uniform_disk/libuniform_disk.dylib"; // module path
	int SUBMODULE = 0; // submodule index (meandr monocycle)

	ModuleManager mng = ModuleManager(NULL);
	mng.load_module(MODULE, R, A0, tau0, EPS, MU);
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;

	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau0] (double vt) {return Function::gauss(vt,tau0);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return linear;
}

vector<SpaceInterval64> arguments (double tau0, bool swap_axis)
{
	vector<SpaceInterval64> data;

	double x = 0;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= 5; z += 0.05) {
			double rho = sqrt(x*x + y*y);
			double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau0 + sqrt((rho+R)*(rho+R) + z*z) + 0.01;
			if (swap_axis) data.emplace_back(y, x, z, from, to);
			else data.emplace_back(x, y, z, from, to);
		}
	}

	return data;
}

void plot (const vector<vector<double>>& data, double tau0, bool swap_axis)
{
	string script = "W_" + str_of(tau0) + (swap_axis?"_xz" :"_yz") + ".gnp";
	GnuPlot* plot = new GnuPlot(script);
	plot->set_colormap(Colormap::parula);
    plot->plot_colormap(data, !swap_axis, 2);
	delete plot;
}

void plot_energy_distribution (double tau0, bool swap_axis)
{
	AbstractField* model = create_model (tau0);
	vector<SpaceInterval64> data = arguments(tau0, swap_axis);

	Cluster cluster(6);
	std::vector<double> res(data.size());
	cluster.start(data, res, [model] (const SpaceInterval64 &entry) {
		return model->energy_cart(entry.q1, entry.q2, entry.q3, entry.from, entry.to);
	});

	while (cluster.progress() != data.size()) {
		double persent = static_cast<double>(cluster.progress()) / data.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	cout << endl;

	// std::cout << "Wainting...\n";
	// cluster.wait();

    vector<vector<double>> ans;
    ans.reserve(data.size());
    for (auto i = 0u; i < data.size(); i++) {
        vector<double> v{ data[i].q1, data[i].q2, data[i].q3, data[i].from, data[i].to, res[i] };
        ans.emplace_back(move(v));
    }

	plot(ans, tau0, swap_axis);
	delete model;
}

int main ()
{
    plot_energy_distribution(2.0, false); // YZ
	plot_energy_distribution(2.0, true ); // XZ
    return 0;
}
