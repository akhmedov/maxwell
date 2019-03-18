//
//  snr_dataset.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "dataset.hpp"

#include <vector>
#include <iostream>
#include <functional>
using namespace std;

static const int 	RADIX  = 3;
static const float 	DCYCLE = 0.5;
static const float 	NPOWER = 10;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity
static const double TAU0 = 1; // duration of signals

static vector<AbstractField*> SIGNAL;

AbstractField* arbitrary_signal (const function<double(double)>& shape) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 1; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, NAN, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	LinearCurrent* source = mng.get_module(mng.get_loaded()[SUBMODULE]).source;
	LinearMedium* medium = mng.get_module(mng.get_loaded()[SUBMODULE]).medium;
	AbstractField* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;

	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth(shape);
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, (LinearField*) linear, NULL);
	return duhamel;
}

void append (int id, Dataset* ds)
{
    vector<double> time, func;

    double rho = 0, phi = 0, z = 2;
    vector<double> space {rho,phi,z};

    float from = (rho > R) ? sqrt((rho-R)*(rho-R) + z*z) : z;
	if (from - 0.01 > 0) from -= 0.01;
	float to = TAU0 + sqrt((rho+R)*(rho+R) + z*z) + 0.01;

    for (float ct = from; ct < to; ct += 0.05) {
        time.push_back(ct);
        func.push_back(SIGNAL.at(id)->electric_x(ct,rho,phi,z));
    }

    ds->append(id,space,time,func);
}

int main ()
{
    vector<function<double(double)>> SHAPE = {
        [] (double time) {return Function::gauss_perp(time,TAU0,1);},
        [] (double time) {return Function::sinc(time,TAU0);}
    };

	Dataset* ds = new Dataset(RADIX, DCYCLE, NPOWER);
    for (auto& i : SHAPE) SIGNAL.push_back(arbitrary_signal(i));
    for (int spark = 0; spark < 10; spark++) append(spark % 2 ? 1 : 0, ds);
    Dataset::serialize("one-by-one.json", ds->get_dataset("one-by-one"), false);
    return 0;
}
