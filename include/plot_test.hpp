//
//  plot_test.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 31.01.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef plot_test_hpp
#define plot_test_hpp

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#include "config.hpp"
#include "manager.hpp"
#include "function.hpp"
#include "gnu_plot.hpp"
#include "integral.hpp"
#include "gnu_plot.hpp"
#include "phys_math.hpp"
#include "kerr_amendment.hpp"
#include "linear_duhamel.hpp"
#include "uniform_disk_current.hpp"

struct PlotTest : private Math, private Integral, private KerrAmendment {

	static void signal_spectr ();
	static void plot_energy (double z);
	static void Ex_derivative (double tau);
	static void plot_source ();
	static void kerr_under_integral_from_nu (const std::vector<double> &r, const std::vector<double> &r_perp, double R = 1);
	static void nonlinear_current (double phi = 0, double z = 0);
	static void I1_numeric_integral (std::size_t points = 1e4, double limits = 1e5);
	static void numeric_perp (double omega = 1);
	static void I1_I2_versus (double rho, double z, double R = 1);
	static void I1_time_partder (double rho, double z, double R = 1);
	static void I2_time_partder (double rho, double z, double R = 1);
	static void kerr_ammend_undeintegral (double ct_perp, double phi_perp, double rho_perp, double z_perp, double R = 1);
	static void set_options ();

private:
	static Config* global_conf;
};

#endif /* plot_test_hpp */