//
//  plot_mangnitude_surface.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"

using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = R; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

AbstractField<Point::Cylindrical>* load_model () 
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
	AbstractField<Point::Cylindrical>* model = load_model();

	double x_from = 0, x_to = 1.5, y = 0, z = 2;
	double from = model->observed_from(Point::Cylindrical::convert(Point::Cartesian3D(x_from,y,z)));
	double to   = model->observed_to(Point::Cylindrical::convert(Point::Cartesian3D(x_to,y,z)));

    vector<Point::SpaceTime<Point::Cylindrical>> args;
    for (double x = x_from; x < x_to; x += 0.05) {
        for (double ct = from; ct < to; ct += 0.05) {
			Point::Cylindrical point = Point::Cylindrical::convert(Point::Cartesian3D(x,y,z));
			args.emplace_back(ct, point);
		}
	}

	vector<double> Z(args.size());
	CalculationManager cluster(6);
	cluster.start(args, Z, model, &AbstractField<Point::Cylindrical>::electric_x);
	cluster.wait();

	vector<pair<double,double>> XY;
	for (size_t i = 0; i < args.size(); i++) {
		Point::Cartesian3D pt = Point::Cartesian3D::convert(args[i]);
        XY.emplace_back(args[i].ct(), pt.x());
    }

	PyPlotManager plot = PyPlotManager("Ex_ct_x.py");
	plot.set_title("Ex(ct,x)");
	plot.set_ox_label("ct, R");
	plot.set_oy_label("x, R");
	plot.set_colormap(ScriptManager::Colormap::jet); // grey, jet, hot, coolwarm
	plot.plot3d(XY, Z);
    return 0;
}
