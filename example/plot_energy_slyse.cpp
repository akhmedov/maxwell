//
//  plot_energy_slyse.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"
#include "module_manager.hpp"

#include <vector>
#include <string> // string() to_string()
#include <utility> // swap()
#include <iostream> // cout
using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

auto str_of = [] (double val) { return to_string(val).substr(0,4); };

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
	free_shape->set_time_depth([tau0] (double vt) {return Function::gauss_perp(vt,tau0,1);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return duhamel;
}

vector<vector<double>> arguments (double tau0, double z)
{
	double range = z;
	vector<vector<double>> data;

	for (double x = -range; x <= range; x += 0.05) {
		for (double y = -range; y <= range; y += 0.05) {
			double rho = sqrt(x*x + y*y);
			double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau0 + sqrt((rho+R)*(rho+R) + z*z) + 0.01;
			vector<double> tmp { x, y, z, from, to };
			data.emplace_back(move(tmp));
		}
	}

	return data;
}

void plot (const vector<vector<double>>& arg, const vector<double>& res, double tau0, double z)
{
    vector<vector<double>> data;
    data.reserve(arg.size());
    for (auto i = 0u; i < arg.size(); i++) {
        vector<double> tmp{ arg[i][0], arg[i][1], arg[i][2], arg[i][3], arg[i][4], res[i] };
        data.emplace_back(move(tmp));
    }

	string script = str_of(tau0) + "_" + str_of(z) + ".gnp" ;
	GnuPlot* plot = new GnuPlot(script);
	plot->set_colormap(Colormap::gray);
	plot->plot_colormap(data, 0, 1);
}

void print_max (const vector<vector<double>>& args, const vector<double>& res)
{
	size_t idx = max_element(res.begin(), res.end()) - res.begin();

	cout << "  x = " << str_of(args[idx][0]);
	cout << "; y = " << str_of(args[idx][1]);
	cout << "; z = " << str_of(args[idx][2]);
	cout << "W(x,y,z) = " << str_of(res[idx]) << endl;
}

void plot_energy_slyse (double tau0, double z)
{
	AbstractField* model = arbitrary_signal(tau0);
	auto args = arguments (tau0, z);
	vector<double> res(args.size());

	CalculationManager cluster(6);
	cluster.start(args, res, [model] (const auto& arg) {
		return model->energy_cart(arg[0], arg[1], arg[2], arg[3], arg[4]);
	});

	// progress of calculation in persent
	while (cluster.progress() != args.size()) {
		double persent = static_cast<double>(cluster.progress()) / args.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	// std::cout << "Wainting...\n";
	// cluster.wait();

	print_max(args, res);
	plot(args, res, tau0, z);
}

int main ()
{
    plot_energy_slyse(0,2);
    return 0;
}
