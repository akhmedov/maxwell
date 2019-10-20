//
//  plot_energy_distribution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

#include <cstdlib> // atof
#include <vector>
#include <typeinfo> // typeid
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

int main (int argc, char* argv[])
{
	if (argc != 2) { 
		cout << "Fail: wrong argument" << endl;
		return 0;
	}

	float tau0 = atof(argv[1]);

	AbstractField<Point::Cylindrical>* model = rectangular_shape(tau0);

	vector<Point::Cylindrical> arg;
	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= 5*R; z += 0.05) {
			Point::Cartesian3D cart = Point::Cartesian3D(0,y,z);
			Point::Cylindrical cylindr = Point::Cylindrical::convert(cart);
			arg.push_back(cylindr);
		}
	}

	vector<double> res(arg.size());
	CalculationManager cluster(4);
	cluster.start(arg, res, model, &AbstractField<Point::Cylindrical>::energy_e);

	while (cluster.progress() != arg.size()) {
		double persent = static_cast<double>(cluster.progress()) / arg.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;
	delete model;

	vector<pair<double,double>> plot_arg;
	plot_arg.reserve(arg.size());
	for (auto i = 0u; i < arg.size(); i++) {
		Point::Cartesian3D point = Point::Cartesian3D::convert(arg[i]);
		plot_arg.emplace_back(point.y(), point.z());
	}

	PyPlotManager plot = PyPlotManager("We_rect_"+ to_string(tau0) +"_YZ.py");
	plot.set_title("Rect pulse duration: " + to_string(tau0));
	plot.set_ox_label("OY, R");
	plot.set_oy_label("OZ, R");
	plot.set_colormap(ScriptManager::Colormap::jet); // grey, jet, hot, coolwarm
	plot.colormap(plot_arg, res);
    return 0;
}
