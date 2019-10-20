//
//  emp_shape_plot.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 09.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

#include <vector>
#include <iostream>
using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = 0.4 * R; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

AbstractField<Point::Cylindrical>* create_model () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 2; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, TAU0, EPS, MU);

	if (!loaded) throw std::logic_error("Library loading error");
	AbstractField<Point::Cylindrical>* tr = mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	auto shape = [] (double ct) { return Function::gauss_perp(ct,TAU0,1); };
	return new DuhamelSuperpose<CylindricalField,Point::Cylindrical>(tr, TAU0, shape, NULL);
}

int main ()
{
    AbstractField<Point::Cylindrical>* model = create_model();

	Point::Cartesian3D cart = Point::Cartesian3D(0.5,0,2);
	Point::Cylindrical cyln = Point::Cylindrical::convert(cart);
	double from = model->observed_from(cyln) - 0.1;
	double to = model->observed_to(cyln) + 0.1;
	Point::SpaceTime<Point::Cylindrical> event {cyln};

	vector<double> arg, fnc;
	for (event.ct() = 1.8; event.ct() <= 2.8; event.ct() += 0.01) {
		fnc.push_back(model->electric_x(event));
		arg.push_back(event.ct());
	}

	PyPlotManager plot = PyPlotManager("emp_shape.py");
	plot.set_ox_label("ct, R");
	plot.set_oy_label("Ex, A/m");
	plot.set_colormap(ScriptManager::Colormap::grey);
	plot.plot2d(arg, fnc);

	delete model;
    return 0;
}
