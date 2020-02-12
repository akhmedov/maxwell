//
//  nonlinear_current.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = R; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity
static const double CHI3 = 1; // relative dielectric pirmativity

AbstractField<Point::Cylindrical>* create_model () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (UniformDisk.NonlinearKerrAmendment)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, TAU0, EPS, MU);

	if (!loaded) throw std::logic_error("Library loading error");
	AbstractField<Point::Cylindrical>* kerr = mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return kerr;
}


vector<double> nonlinear_current (AbstractField<Point::Cylindrical>* model, const vector<Point::ModalSpaceTime<Point::Cylindrical>>& events)
{
	vector<double> res(events.size());
	CalculationManager cluster(5);
	cluster.start(events, res, model, &AbstractField<Point::Cylindrical>::modal_jm);

	auto str_of = [] (double val) { return to_string(val).substr(0,4); };
	while (cluster.progress() != events.size()) {
		double persent = static_cast<double>(cluster.progress()) / events.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	// cluster.wait();
	delete model;
	return res;
}

void plot_colormap (const vector<Point::ModalSpaceTime<Point::Cylindrical>>& events, const vector<double>& res)
{
	vector<pair<double,double>> args;
	args.reserve(events.size());
	for (auto i = 0u; i < events.size(); i++) {
		double ct_z = events[i].ct() - events[i].z();
		double nu = events[i].nu();
		args.emplace_back(nu, ct_z);
	}

	PyPlotManager plot = PyPlotManager("kerr_source.py");
	plot.set_title("jm, V/m");
	plot.set_ox_label("nu");
	plot.set_oy_label("ct_z, R");
	plot.set_colormap(ScriptManager::Colormap::coolwarm); // grey, jet, hot, coolwarm
	plot.colormap(args, res);
}

int main ()
{
	AbstractField<Point::Cylindrical>* kerr = create_model();
	vector<Point::ModalSpaceTime<Point::Cylindrical>> events;

	for (double ct = 2; ct < 2.3; ct += 0.01) {
		for (double nu = 0; nu < 20; nu += 0.01) {
			Point::ModalSpaceTime<Point::Cylindrical> observer = {1, nu, ct, 0.0, 0.0, 2.0};
			events.push_back(observer);
		}
	}

	vector<double> res = nonlinear_current(kerr, events);
	plot_colormap(events, res);

    return 0;
}
