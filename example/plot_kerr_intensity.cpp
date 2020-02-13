//
//  plot_kerr_intensity.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 15.02.20.
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

vector<double> calculate_evolution (AbstractField<Point::Cylindrical>* model, const vector<Point::ModalSpaceTime<Point::Cylindrical>>& events)
{
	vector<double> res(events.size());
	DatabaseCalculationManager cluster(12, "plot_kerr_evolution.example.maxwell");
	cluster.start(events, res, model, &AbstractField<Point::Cylindrical>::modal_vmh);

	cluster.wait();

	/* auto str_of = [] (double val) { return to_string(val).substr(0,4); };
	while (cluster.progress() != events.size()) {
		double persent = static_cast<double>(cluster.progress()) / events.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl; */

	return res;
}

void plot_intensity (const vector<Point::SpaceTime<Point::Cylindrical>>& event, const vector<double>& field)
{
	vector<double> time;
	for (auto i : event) time.push_back(i.ct());

	PyPlotManager plot = PyPlotManager("kerr_intensity.py");
	plot.set_ox_label("ct, R");
	plot.set_oy_label("Ex, V/m");
	plot.set_colormap(ScriptManager::Colormap::grey);
	plot.plot2d(time, field);
}

int main ()
{
	Point::Cylindrical position(0.0, 0.0, 2.0);

	vector<Point::SpaceTime<Point::Cylindrical>> event;
	for (double ct = 2.00; ct <= 2.60; ct += 0.05) {
		event.emplace_back(ct, position);
	}

	vector<Point::ModalSpaceTime<Point::Cylindrical>> modal;
	for (auto i : event) {
		for (double nu = 0; nu < 69; nu += 0.25) {
			modal.emplace_back(1, nu, i);
			modal.emplace_back(3, nu, i);
		}
	}

	cout << "[info] Argument ready, number of items: " << modal.size() << endl;

	AbstractField<Point::Cylindrical>* kerr = create_model();
	vector<double> res = calculate_evolution(kerr, modal);

	cout << "[info] Evolution coefficients are ready." << endl;

	vector<double> field;
	for (auto i : event) field.push_back(kerr->electric_x(i));

	cout << "[info] Field intensity ready." << endl;

	plot_intensity(event, field);

    return 0;
}
