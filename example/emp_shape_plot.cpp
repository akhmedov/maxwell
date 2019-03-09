//
//  emp_shape_plot.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 09.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

using namespace std;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = 2; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

AbstractField* create_model () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, TAU0, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
    cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	    

	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([] (double vt) {return Function::gauss_perp(vt,TAU0,1);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return duhamel;
}

vector<vector<double>> plot_data (AbstractField* model, double x, double y, double z)
{
    double rho = sqrt(x*x + y*y);
    double phi = atan2(y,x);
    double from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : abs(z);
    double to = TAU0 + sqrt((rho+R)*(rho+R) + z*z);

    vector<vector<double>> data;
	for (double ct = from; ct <= to; ct += 0.05) {
        double Ex = model->electric_x(ct, rho, phi, z);
		data.push_back({ct, Ex});
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
    AbstractField* model = create_model();
    auto data = plot_data(model, 1, 1, 1.6);
    plot_script(data, "emp_shape.gnp");

    return 0;
}
