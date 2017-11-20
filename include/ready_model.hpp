//
//  ready_model.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 04.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "uniform_disk_current.hpp"
#include "gnu_plot.hpp"
#include "manager.hpp"
#include "integral.hpp"

#include <vector>
#include <tuple> // std::make_tuple
#include <limits> // std::numeric_limits<size_t>::max

struct ReadyModel {
	/* 01 */ static void linear_Hy_from_ct (double rho = 0, double phi = 0, double z = 0.01);
	/* 02 */ static void linear_Ex_from_ct (double rho = 0, double phi = 0, double z = 0);
	/* 03 */ static void linear_Hy_from_z (double ct = 100, double rho = 0, double phi = 0);
	/* 04 */ static void linear_Hy_from_ct_z (double rho = 0, double phi = 0);
	/* 05 */ static void linear_Ex_from_ct_z (double rho = 0, double phi = 0);
	/* 06 */ static void linear_Ex_from_rho_phi (double ct = 0.5, double z = 0);
	/* 07 */ static void linear_Ex_from_ct_rho (double phi = 0, double z = 2);
	/* 08 */ static void missile_effect_length ();
	/* 09 */ static void magnetic_static_magnitude (double eps = 10e-10, double rho = 0, double phi = 0);
	/* 10 */ static void magnetic_static_terms (double rho = 0, double phi = 0, double z = 0.9);
	/* 11 */ static void inegration_compare_i1 ();
	/* 12 */ static void inegration_compare_i2 ();
};
