//
//  plot_kerr_evolution.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 09.02.20.
//  Copyright © 2020 Rolan Akhmedov. All rights reserved.
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

vector<double> nonlinear_evolution (AbstractField<Point::Cylindrical>* model, const vector<Point::ModalSpaceTime<Point::Cylindrical>>& events)
{
	vector<double> res(events.size());
	CalculationManager cluster(4);
	cluster.start(events, res, model, &AbstractField<Point::Cylindrical>::modal_vmh);

	auto str_of = [] (double val) { return to_string(val).substr(0,4); };
	while (cluster.progress() != events.size()) {
		double persent = static_cast<double>(cluster.progress()) / events.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	delete model;
	return res;
}

void plot_evolution (const vector<Point::ModalSpaceTime<Point::Cylindrical>>& events, const vector<double>& res)
{
	vector<double> arg;
	arg.reserve(events.size());
	for (auto i = 0u; i < events.size(); i++)
		arg.push_back(events[i].nu());

	PyPlotManager plot = PyPlotManager("kerr_evolution.py");
	plot.set_ox_label("nu");
	plot.set_oy_label("Vmh");
	plot.set_colormap(ScriptManager::Colormap::grey); // grey, jet, hot, coolwarm
	plot.plot2d(arg, res);
}

int main ()
{
	AbstractField<Point::Cylindrical>* kerr = create_model();
	vector<Point::ModalSpaceTime<Point::Cylindrical>> events;

	for (double nu = 0; nu < 10; nu += 1) {
		Point::ModalSpaceTime<Point::Cylindrical> observer = {1, nu, 2.1, 0.0, 0.0, 2.0};
		events.push_back(observer);
	}

	cout << "Argument ready, number of items: " << events.size() << endl;

	vector<double> res = nonlinear_evolution(kerr, events);
	plot_evolution(events, res);

    return 0;
}
