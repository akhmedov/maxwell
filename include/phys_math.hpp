//
//  phys_math.hpp
//  Evolution
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef phys_math_hpp
#define phys_math_hpp

#include <gmp.h>
#include <gmpxx.h>

#include <string>
#include <cmath>
#include <cfloat>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <functional>

#define DERIVATIVE_STEP 1e-10;

typedef std::pair<std::size_t,std::size_t> Binom;

struct Math {
	static double deg2rad (double deg_angle);
	static bool compare (double val1, double val2);
	static double kronecker_delta (double arg, double param);
	static double heaviside_sfunc (double arg, double param, bool mirror = false);
	static mpf_class yacobi_polinom (std::size_t term, mpf_class arg, std::size_t alpha = 1, std::size_t beta = 0);
	static double yacobi_polinom_recur (std::size_t term, double arg, std::size_t alpha = 1, std::size_t beta = 0);
	static double yacobi_polinom_init (std::size_t term, double arg, std::size_t alpha = 1, std::size_t beta = 0);
	static std::size_t two_in_pow (std::size_t exponent);
	static std::size_t binomial (std::size_t n, std::size_t k);
	static std::size_t product (std::size_t first, std::size_t second);
	static double divide (std::size_t one, std::size_t two);
	static float inv_sqrt (float arg);
	static double derivat3 ( std::function<double(double)>, double arg );
	static double derivat4 ( std::function<double(double)>, double arg );
	static double lommel (std::size_t terms, int order, double W, double Z);
protected:
	static std::size_t next_prime(std::size_t prime, std::size_t search_limit);
	static std::size_t binom_prod (std::size_t n, std::size_t m);
	static void divide3to1 (std::size_t &x, std::size_t &y, std::size_t &z, std::size_t k);
	static const std::size_t MAX_DIVIDE_PRECISION = __SIZE_MAX__ / 10e2;
};

#endif /* phys_math_hpp */
