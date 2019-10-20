//
//  plot_energy_slyse.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

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

vector<Point::Cylindrical> arguments (double z)
{
	double range = z;
	vector<Point::Cylindrical> data;

	for (double x = -range; x <= range; x += 0.05) {
		for (double y = -range; y <= range; y += 0.05) {
			Point::Cartesian3D cart = Point::Cartesian3D(x,y,z);
			Point::Cylindrical cyl = Point::Cylindrical::convert(cart);
			data.push_back(cyl);
		}
	}

	return data;
}

void plot (const vector<Point::Cylindrical>& arg, const vector<double>& res, double tau0, double z)
{
	vector<pair<double,double>> XY;
    XY.reserve(arg.size());
    for (auto i = 0u; i < arg.size(); i++) {
		Point::Cartesian3D point = Point::Cartesian3D::convert(arg[i]);
        vector<double> tmp{ point.x(), point.y(), point.z(), res[i] };
		XY.emplace_back(point.x(), point.y());
    }

	string script = "We_" + str_of(tau0) + "_" + str_of(z) + ".py" ;
	PyPlotManager plot = PyPlotManager(script);
	plot.set_title("Pulse duration: " + to_string(tau0));
	plot.set_ox_label("OX, R");
	plot.set_oy_label("OY, R");
	plot.set_colormap(ScriptManager::Colormap::coolwarm); // grey, jet, hot, coolwarm
	plot.colormap(XY, res);
}

void print_max (const vector<Point::Cylindrical>& args, const vector<double>& res)
{
	auto idx = max_element(res.begin(), res.end()) - res.begin();
	cout << "x,y,z = " << args[idx].to_str() << endl;
	cout << "W(x,y,z) = " << std::to_string(res[idx]) << endl;
}

void plot_energy_slyse (double tau0, double z)
{
	AbstractField<Point::Cylindrical>* model = arbitrary_signal(tau0);
	vector<Point::Cylindrical> args = arguments (z);
	vector<double> res(args.size());

	CalculationManager cluster(6);
	cluster.start(args, res, model, &AbstractField<Point::Cylindrical>::energy_e);

	// progress of calculation in persent
	while (cluster.progress() != args.size()) {
		double persent = static_cast<double>(cluster.progress()) / args.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	print_max(args, res);
	plot(args, res, tau0, z);
}

int main ()
{
    plot_energy_slyse(0,2);
    return 0;
}
