//
//  plot_mangnitude_surface.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

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

vector<vector<double>> plot_data (AbstractField<Point::Cylindrical>* model, double x_from, double x_to, double y, double z)
{
	double from = model->observed_from(Point::Cylindrical::convert(Point::Cartesian3D(x_from,y,z)));
	double to   = model->observed_to(Point::Cylindrical::convert(Point::Cartesian3D(x_to,y,z)));

    vector<Point::SpaceTime<Point::Cylindrical>> args;
    for (double x = x_from; x < x_to; x += 0.05) {
        for (double ct = from; ct < to; ct += 0.05) {
			Point::Cylindrical point = Point::Cylindrical::convert(Point::Cartesian3D(x,y,z));
			args.emplace_back(ct, point);
		}
	}

	vector<double> res(args.size());
	CalculationManager cluster(6);
	cluster.start(args, res, model, &AbstractField<Point::Cylindrical>::electric_x);
	cluster.wait();

	vector<vector<double>> data;
	for (size_t i = 0; i < args.size(); i++) {
		Point::Cartesian3D pt = Point::Cartesian3D::convert(args[i]);
        data.push_back({ res[i], args[i].ct(), pt.x() });
    }

	return data;
}

void plot (const vector<vector<double>>& data)
{
	GnuPlot plot = GnuPlot("emp_surface.gnp");
    plot.plot3d(data);
}

int main ()
{
	AbstractField<Point::Cylindrical>* model = load_model();
	vector<vector<double>> data = plot_data(model, 0, 1.5, 0, 2);
	plot(data);
    return 0;
}
