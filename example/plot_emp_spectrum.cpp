//
//  plot_emp_spectrum.hpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 11.10.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "phys_math.hpp"
#include "pyplot_manager.hpp"

#include <vector>
#include <iostream>
using namespace std;

static const double time_resolution = 0.01;

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

	Point::Cartesian3D cart = Point::Cartesian3D(0.5,0.5,2);
	Point::Cylindrical cyln = Point::Cylindrical::convert(cart);
	Point::SpaceTime<Point::Cylindrical> event {cyln};

	vector<double> arg;
	vector<complex<double>> fnc;
	for (event.ct() = model->observed_from(cyln); event.ct() <= model->observed_to(cyln); event.ct() += time_resolution) {
		fnc.emplace_back(model->electric_x(event), 0);
		arg.push_back(event.ct());
	}

    fnc = Math::fft(fnc);
    for_each(fnc.begin(), fnc.end(), [](complex<double> &n){ n *= time_resolution; });

    vector<double> spectr(fnc.size());
    transform(fnc.begin(), fnc.end(), spectr.begin(), [] (complex<double> z) { return z.real(); });
    arg = Math::fft_arg(fnc.size(), time_resolution);

	fnc.push_back(*fnc.begin());
	arg.push_back(-*arg.begin());

	vector<double> fnc_re, fnc_im;
	fnc_re.reserve(fnc.size());
	fnc_im.reserve(fnc.size());
	for (auto i : fnc) {
		fnc_re.push_back(i.real());
		fnc_im.push_back(i.imag());
	}

	PyPlotManager plot = PyPlotManager("emp_spectrum.py");
	plot.set_ox_label("f/R, GHz");
	plot.set_oy_label("????");
	plot.set_colormap(ScriptManager::Colormap::grey);
	plot.plot2d({arg, arg}, {fnc_re, fnc_im}, {"Re(fz)", "Im(fz)"});

	delete model;
    return 0;
}
