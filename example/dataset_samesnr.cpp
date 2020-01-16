//
//  snr_dataset.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "dataset.hpp"

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <functional>
using namespace std;

static const int 	RADIX  = 4;
static const float 	DCYCLE = 150;
static const float 	NPOWER = 10;

static const double R  = 1; // disk radius
static const double A0 = 1; // max current magnitude
static const double TAU0 = R; // duration of signals

static const double MU  = 1; // relative magnetic permatiity
static const double EPS = 1; // relative dielectric pirmativity

static vector<AbstractField<Point::Cylindrical>*> SIGNAL;

AbstractField<Point::Cylindrical>* arbitrary_signal (const function<double(double)>& shape) 
{
	string MODULE_PATH = "module/uniform_disk"; // module dir path
	string MODULE_NAME = "uniform_disk"; // library name
	int SUBMODULE = 2; // submodule index (trancient responce)

	ModuleManager mng = ModuleManager(NULL);
	bool loaded = mng.load_module(MODULE_PATH, MODULE_NAME, R, A0, NAN, EPS, MU);
	if (!loaded) throw std::logic_error("Library loading error");
	AbstractField<Point::Cylindrical>* linear = mng.get_module(mng.get_loaded()[SUBMODULE]).field_cyl_arg;
	cout << "Submodule loaded: " << mng.get_loaded()[SUBMODULE] << endl;
	return new DuhamelSuperpose<CylindricalField,Point::Cylindrical>(linear, TAU0, shape);
}

void append (double rho, double phi, double z, int id, Dataset::Dataset* ds)
{
    vector<double> time, func;

    std::vector<double> coord {rho,phi,z};
    Point::Cylindrical observer {rho,phi,z};

    for (float ct = SIGNAL.at(id)->observed_from(observer); ct < SIGNAL.at(id)->observed_to(observer); ct += 0.02) {
        Point::SpaceTime<Point::Cylindrical> event{ct,rho,phi,z};
        func.push_back(SIGNAL.at(id)->electric_x(event));
    }

    Dataset::Annotation annotation = Dataset::Annotation(id, 0.02, coord);
    ds->append(func, annotation);
}

int main ()
{
    vector<function<double(double)>> SHAPE = {
        [] (double time) {return Function::gauss(time,TAU0);},
        [] (double time) {return 0.5 * Function::gauss_perp(time,TAU0,1);},
        [] (double time) {return Function::sinc(time,TAU0);}
    };

    SIGNAL.push_back(new ZeroField<Point::Cylindrical>());
    for (auto& i : SHAPE) SIGNAL.push_back(arbitrary_signal(i));

    std::random_device random;
    std::mt19937 generator(random());
    std::uniform_int_distribution<std::size_t> signal_id_distr(0, SIGNAL.size() - 1);
    std::uniform_real_distribution<double> rho_distr(0, R);
    std::uniform_real_distribution<double> phi_distr(0, M_PI_2);
    std::uniform_real_distribution<double> z_distr(R, 3*R);

    vector<string> class_label = {"void", "gauss", "gauss_perp", "sinc"};
	Dataset::Dataset* dataset = new Dataset::Dataset(RADIX, DCYCLE, NPOWER, class_label);
    for (int spark = 0; spark < 10000; spark++) {
        if (spark % 100 == 0) cout << "Sparks ready:" << spark << endl;
        size_t signal_id = signal_id_distr(generator);
        double rho = rho_distr(generator);
        double phi = phi_distr(generator);
        double z = z_distr(generator);
        append(rho, phi, z, signal_id, dataset);
    }

    dataset->serialize_to_json("one-by-one.json", false);
    return 0;
}
