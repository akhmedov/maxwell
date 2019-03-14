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
static const double TAU0 = 2; // duration of exitation

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

auto str_of = [] (double val) { return to_string(val).substr(0,4); };

AbstractField* create_model () 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, TAU0, EPS, MU);
	if (!loaded) throw logic_error("Library loading error");
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
    cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	    
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([] (double vt) {return Function::gauss_perp(vt,TAU0,1);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return duhamel;
}

vector<vector<double>> argument (double x_from, double x_to, double y, double z)
{
	double min_rho = sqrt(x_from * x_from + y*y);
	double max_rho = sqrt(x_to * x_to + y*y);

	double from = (min_rho > R) ? sqrt((min_rho-R)*(min_rho-R) + z*z) : abs(z);
	double to = TAU0 + sqrt((max_rho+R)*(max_rho+R) + z*z);

    vector<vector<double>> data;

    for (double x = x_from; x < x_to; x += 0.05) {

		double rho = sqrt(x*x + y*y);
		double phi = atan2(y,x);

        for (double ct = from; ct < to; ct += 0.05)
			data.push_back({ct,rho,phi,z});
    }

	return data;
}

vector<vector<double>> calculator (AbstractField* model, const vector<vector<double>>& args)
{
	vector<double> res(args.size());

	CalculationManager cluster(6);
	cluster.start(args, res, [model] (const auto& item) {
		return model->electric_x(item[0], item[1], item[2], item[3]);
	});

	// progress of calculation in persent
	while (cluster.progress() != args.size()) {
		double persent = static_cast<double>(cluster.progress()) / args.size();
		cout << '\r' << "Evaluation progress: " << str_of(100 * persent) << '%';
		cout.flush();
		this_thread::sleep_for(chrono::milliseconds(100));
	}
	cout << endl;

	vector<vector<double>> data;
	for (size_t i = 0; i < args.size(); i++) {
		double x  = args[i][1] * cos(args[i][2]);
		// double y  = args[i][1] * sin(args[i][2]);
        data.push_back({ res[i], args[i][0], x });
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
	AbstractField* model = create_model();
	auto args = argument(0, 1.5, 0, 2);
	auto data = calculator (model, args);
	plot(data);
    return 0;
}