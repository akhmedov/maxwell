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

AbstractField<Point::Cylindrical>* trancient_recponce () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 2; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, 0.0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
}

AbstractField<Point::Cylindrical>* rectangular_shape (double tau0) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 0; // submodule index (meandr monocycle)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, tau0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
}

AbstractField<Point::Cylindrical>* arbitrary_signal (double tau0) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 2; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, tau0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	AbstractField<Point::Cylindrical>* tr = mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	auto shape = [tau0] (double ct) { return Function::gauss_perp(ct,tau0,1); };
	return new DuhamelSuperpose<CylindricalField,Point::Cylindrical>(tr, tau0, shape, NULL);
}

vector<Point::Cylindrical> arguments (bool swap_axis)
{
	Point::Cartesian3D cart;
	vector<Point::Cylindrical> data;

	double x = 0;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= 5*R; z += 0.05) {
			if (swap_axis) cart = Point::Cartesian3D(y,x,z);
			else cart = Point::Cartesian3D(x,y,z);
			data.push_back(Point::Cylindrical::convert(cart));
		}
	}

	return data;
}

void plot (const vector<Point::Cylindrical>& arg, const vector<double>& res, double tau0, bool swap_axis)
{
    vector<vector<double>> data;
    data.reserve(arg.size());
    for (auto i = 0u; i < arg.size(); i++) {
		Point::Cartesian3D point = Point::Cartesian3D::convert(arg[i]);
        vector<double> tmp{ point.x(), point.y(), point.z(), res[i] };
        data.emplace_back(move(tmp));
    }

	string script = "W_" + str_of(tau0) + (swap_axis ? "_xz" : "_yz") + ".gnp";
	GnuPlot plot = GnuPlot(script);
	plot.set_colormap(Colormap::gray);
    plot.plot_colormap(data, !swap_axis, 2);
}

void plot_energy_distribution (double tau0, bool swap_axis)
{
	AbstractField<Point::Cylindrical>* model = trancient_recponce();
	vector<Point::Cylindrical> arg = arguments(swap_axis);
	vector<double> res(arg.size());

	CalculationManager cluster(4);
	cluster.start(arg, res, model, &AbstractField<Point::Cylindrical>::energy_e);

	// progress of calculation in persent
	while (cluster.progress() != arg.size()) {
		double persent = static_cast<double>(cluster.progress()) / arg.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

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
