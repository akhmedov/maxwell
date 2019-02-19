//
//  unit_test.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 08.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef unit_test_hpp
#define unit_test_hpp

#include <gmp.h>
#include <gmpxx.h>

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#include <cfloat> // size_t
#include <iomanip> // std::fixed std::setprecision(4)
#include <iostream> // cout cerr
#include <string> // compare (==)
#include <vector>
#include <fstream>

#include "dataset.hpp"
#include "phys_math.hpp"
#include "integral.hpp"
#include "gnu_plot.hpp"
#include "function.hpp"
#include "linear_duhamel.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

struct UnitTest : private GnuPlot, private Math, private Integral, private KerrAmendment {
	
	static bool hardcode_dataset ();
	static bool random_dataset ();
	static bool snr_dataset ();
	static bool next_prime ();
	static bool binom_prod ();
	static bool bessel_perp ();
	static bool simpson_runge ();
	static bool simpson_I1 ();
	static bool simpson_I2 ();
	static bool simpson2d ();
	static bool simpson3d ();
	static bool invers_sqrt ();
	static bool field_getters ();
	static bool real_convolution ();
	static bool monte_carlo_vector ();
	static bool yacobi_pol_property ();
	static bool monte_carlo_integral ();
	static bool monte_carlo_improper ();
	static bool imptoper_int_bessel ();
	static bool simpson_runge_jump ();
	static bool yacobi_pol_compare ();
	static bool vector_to_matrix ();
	static bool simpson_runge_2d ();
	static bool I1_time_partder ();
	static bool I2_time_partder ();
	static bool simpson2d_line ();
	static bool simpson2d_line2 ();
	static bool list2grid ();
};

#endif /* unit_test_hpp */
