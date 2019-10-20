//
//  nonlinear_current.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "pyplot_manager.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

void nonlinear_current (double phi, double z)
{
	double R = 1;
	double eps_r = 1;
	double mu_r = 1;
	double xi3 = 0;
	// double inv_em_relation = NonlinearMedium::EPS0 * eps_r / NonlinearMedium::MU0 * mu_r;
	double A0 = 1; // * inv_em_relation;

	Homogeneous* linear_medium = new Homogeneous(mu_r, eps_r);
	KerrMedium* kerr_medium = new KerrMedium(mu_r, eps_r, xi3, 0);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, linear_medium);
	KerrAmendment* non_linear = new KerrAmendment(linear, kerr_medium, source);

	vector<pair<double,double>> args;
	vector<double> func;
	for (double rho = 0.0; rho <= 1.5; rho += 0.01) {
		for (double ct = 0.4; ct <= 1.5; ct += 0.01) {
			if (ct - z <= 0) {
				plot_data.push_back({0,ct,rho});
			} else {
				double jx = non_linear->current_x(ct,rho,phi,z);
				args.emplace_back(ct,rho);
				func.push_back(jx);
			}
		}
	}

	PyPlotManager plot = PyPlotManager("kerr_source.py");
	plot.set_title("jx`(ct,rho), V/m");
	plot.set_ox_label("ct, R");
	plot.set_oy_label("rho, R");
	plot.set_colormap(ScriptManager::Colormap::coolwarm); // grey, jet, hot, coolwarm
	plot.plot3d(args, func);
}

int main ()
{
    nonlinear_current(0,2);
    return 0;
}
