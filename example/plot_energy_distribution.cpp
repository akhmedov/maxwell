//
//  plot_energy_distribution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <vector>
#include <string> // string() to_string()
#include <utility> // swap()
#include <iostream> // cout
using namespace std;

double R  = 1; // disk radius
double A0 = 1; // max current magnitude

double MU  = 1; // relative magnetic permatiity
double EPS = 1; // relative dielectric pirmativity

auto str_of = [] (double val) { return to_string(val).substr(0,4); };

AbstractField* trancient_recponce () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, 0.0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return mng.get_module(mng.get_loaded()[SUBMODULE]).field;
}

AbstractField* rectangular_shape (double tau0) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 0; // submodule index (meandr monocycle)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, tau0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return mng.get_module(mng.get_loaded()[SUBMODULE]).field;
}

AbstractField* arbitrary_signal (double tau0) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, tau0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;

	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau0] (double vt) {return Function::gauss(vt,tau0);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return duhamel;
}

vector<vector<double>> arguments (double tau0, bool swap_axis)
{
	vector<vector<double>> data;

	double x = 0;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= 5*R; z += 0.05) {
			double rho = sqrt(x*x + y*y);
			double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau0 + sqrt((rho+R)*(rho+R) + z*z) + 0.01;
			vector<double> tmp {x, y, z, from, to};
			if (swap_axis) swap(tmp[0],tmp[1]);
			data.emplace_back(move(tmp));
		}
	}

	return data;
}

void plot (const vector<vector<double>>& arg, const vector<double>& res, double tau0, bool swap_axis)
{
    vector<vector<double>> data;
    data.reserve(arg.size());
    for (auto i = 0u; i < arg.size(); i++) {
        vector<double> tmp{ arg[i][0], arg[i][1], arg[i][2], arg[i][3], arg[i][4], res[i] };
        data.emplace_back(move(tmp));
    }

	string script = "W_" + str_of(tau0) + (swap_axis?"_xz" :"_yz") + ".gnp";
	GnuPlot plot = GnuPlot(script);
	plot.set_colormap(Colormap::gray);
    plot.plot_colormap(data, !swap_axis, 2);
}

void plot_energy_distribution (double tau0, bool swap_axis)
{
	AbstractField* model = trancient_recponce();
	auto arg = arguments(tau0, swap_axis);
	vector<double> res(arg.size());

	CalculationManager cluster(3);
	cluster.start(arg, res, [model] (const auto& entry) {
		return model->energy_cart(entry[0], entry[1], entry[2], entry[3], entry[4]);
	});

	// progress of calculation in persent
	while (cluster.progress() != arg.size()) {
		double persent = static_cast<double>(cluster.progress()) / arg.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	// std::cout << "Wainting...\n";
	// cluster.wait();

	plot(arg, res, tau0, swap_axis);
	delete model;
}

int main ()
{
    plot_energy_distribution(1.0, false); // YZ
	plot_energy_distribution(1.0, true ); // XZ

	plot_energy_distribution(2.0, false); // YZ
	plot_energy_distribution(2.0, true ); // XZ

	plot_energy_distribution(3.0, false); // YZ
	plot_energy_distribution(3.0, true ); // XZ

    return 0;
}
