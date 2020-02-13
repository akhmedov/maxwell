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
#include <vector>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <functional>

#define DERIVATIVE_STEP 1e-10;

typedef std::pair<std::size_t,std::size_t> Binom;

struct Math {
	static double log (double base, double deg_angle);
	static double deg2rad (double deg_angle);
	static bool compare (double val1, double val2, double eps = DBL_MIN);
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
	static double simpson_integral (const std::vector<double>& fnc, const std::vector<double>& arg);
	static double lommel (std::size_t terms, int order, double W, double Z);
	static std::vector<std::complex<double>> fft (const std::vector<std::complex<double>> &vec); // DFT for any vec.size()
	static std::vector<double> fft_arg (size_t func_samples, double sampling_rate);
	static std::vector<std::complex<double>> inv_fft (const std::vector<std::complex<double>> &vec); // IDFT for any vec.size()

protected:

	// needed for fft and inv_fft member function
	static void transformRadix2 (std::vector<std::complex<double>> &vec); // Cooley-Tukey FFT algorithm for DFT
	static void transformBluestein (std::vector<std::complex<double>> &vec); // Bluestein's chirp z-transform algorithm
	static void convolve (const std::vector<std::complex<double>> &xvec, const std::vector<std::complex<double>> &yvec, std::vector<std::complex<double>> &outvec); // circular convolution of the given complex vectors
	static size_t reverseBits (size_t x, int n);

	// needed for yacobi polinom calculations
	static std::size_t next_prime(std::size_t prime, std::size_t search_limit);
	static std::size_t binom_prod (std::size_t n, std::size_t m);
	static void divide3to1 (std::size_t &x, std::size_t &y, std::size_t &z, std::size_t k);
	static const std::size_t MAX_DIVIDE_PRECISION = __SIZE_MAX__ / 10e2;
};

#endif /* phys_math_hpp */
