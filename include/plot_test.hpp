//
//  plot_test.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 31.01.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef plot_test_hpp
#define plot_test_hpp

#include "config.hpp"
#include "gnu_plot.hpp"
#include "phys_math.hpp"
#include "integral.hpp"
#include "gnu_plot.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

struct PlotTest : private Math, private Integral, private KerrAmendment {

	static void numeric_perp(double omega = 1);
	static void I1_I2_versus (double rho, double z, double R = 1);
	static void I1_time_partder (double rho, double z, double R = 1);
	static void I2_time_partder (double rho, double z, double R = 1);
	static void kerr_ammend_undeintegral (double ct_perp, double phi_perp, double rho_perp, double z_perp, double R = 1);
	static void set_options ();

private:
	static Config* global_conf;
};

#endif /* plot_test_hpp */