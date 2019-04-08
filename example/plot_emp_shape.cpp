//
//  emp_shape_plot.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 09.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <iostream>
#include <vector>

using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = R; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

AbstractField<Point::Cylindrical>* create_model () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, TAU0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	AbstractField<Point::Cylindrical>* tr = mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	auto shape = [] (double ct) { return Function::gauss_perp(ct,TAU0,1); };
	return new DuhamelSuperpose<CylindricalField,Point::Cylindrical>(tr, TAU0, shape, NULL);
}

template <typename System> vector<vector<double>> plot_arguments (AbstractField<Point::Cylindrical>* model, const System& point)
{
	Point::Cylindrical cyln = Point::Cylindrical::convert(point);
	double from = model->observed_from(cyln) - 0.1;
	double to = model->observed_to(cyln) + 0.1;
	Point::SpaceTime<Point::Cylindrical> event {cyln};

    vector<vector<double>> data;
	for (event.ct() = from; event.ct() <= to; event.ct() += 0.01) {
        double Ex = model->electric_x(event);
		data.push_back({event.ct(), Ex});
    }

	return data;
}

void plot_script (vector<vector<double>> data, string name)
{
	GnuPlot plot = GnuPlot(name);
	plot.set_ox_label("ct, m");
	plot.set_oy_label("Ex, V/m");
	plot.grid_on();
	plot.cage_on();
	plot.plot2d(data);
}

int main ()
{
    AbstractField<Point::Cylindrical>* model = create_model();
    auto data = plot_arguments(model, Point::Cartesian3D(0.5,0.5,2));
    plot_script(data, "emp_shape.gnp");

	delete model;
    return 0;
}
