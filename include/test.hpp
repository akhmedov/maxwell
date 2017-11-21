//
//  test.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 08.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef test_hpp
#define test_hpp

#include <gmp.h>
#include <gmpxx.h>

#include <cfloat> // size_t
#include <iostream> // cout cerr
#include <string> // compare (==)
#include <vector>
#include <fstream>

#include "phys_math.hpp"
#include "ready_model.hpp"
#include "integral.hpp"
#include "gnu_plot.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

struct Test : private Math, private KerrAmendment {
	static bool next_prime ();
	static bool binom_prod ();
	static bool yacobi_pol_property ();
	static bool yacobi_pol_compare ();
	static bool bessel_perp ();
	static bool simpson_I1 ();
	static bool simpson_I2 ();
	static bool invers_sqrt ();
	static bool field_getters ();
	static bool real_convolution ();
};

#endif /* test_hpp */
