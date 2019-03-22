//
//  phys_math.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"

double  Math::log (double base, double arg)
{
	return std::log(arg) / std::log(base);
}
 
double Math::deg2rad (double deg_angle)
{
	return M_PI * deg_angle / 180.0;
}

bool Math::compare(double val1, double val2, double EPS)
{
	double diff = std::abs(val1 - val2);
	if (diff > EPS) return false;
	else return true;
}

double Math::kronecker_delta (double arg, double param)
{
	bool equal = Math::compare(param, arg);
	if (equal) return 1;
	else return 0;
}

// heaviside step function = H(arg - param)
double Math::heaviside_sfunc (double arg, double param, bool mirror)
{
	bool equal = Math::compare(param, arg);
	if (equal) return 0.5;
	if (mirror) {
		if (arg < param) return 1;
		else return 0;
	} else {
		if (arg > param) return 1;
		else return 0;
	}
}

std::size_t Math::two_in_pow (std::size_t exponent)
{
	std::size_t result = 1;
	
	for (std::size_t power = 0; power < exponent; power++)
		result *= 2;
	
	return result;
}

std::size_t Math::binomial (std::size_t n, std::size_t k)
{
	std::size_t res = 1;
	if (k > n - k) k = n - k;
	for (std::size_t i = 1; i <= k; i++, n--) {
		if (res / i > __SIZE_MAX__ / n)
			throw std::overflow_error("Binomial is over of size_t");
		res = res / i * n + res % i * n / i;
	}
	return res;
}

std::size_t Math::product (std::size_t first, std::size_t second)
{
	std::size_t result = first * second;
	bool warning = (result < std::max<std::size_t>(first,second));
	bool is_zero = (first == 0) || (second == 0);
	bool is_one = (first == 1) || (second == 1);
	if (warning && !is_one && !is_zero) 
		throw std::invalid_argument("Product extends size_t value");
	return result;
}

double Math::divide (std::size_t num, std::size_t den)
{
	double result = 0;
	
	std::size_t iter = 1;
	std::size_t devision;
	std::size_t remainder = 42;
	
	while ( (iter <= Math::MAX_DIVIDE_PRECISION) && remainder ) {
		devision = num / den;
		remainder = num % den;
		result += static_cast<double>(devision) / static_cast<double>(iter);
		num = remainder * 10;
		iter *= 10;
	}
	
	return result;
}

mpf_class Math::yacobi_polinom (std::size_t term, mpf_class arg, std::size_t alpha, std::size_t beta)
{
	mpf_class res = 0;

	mpz_class binom1, binom2;
	mpf_class arg_minus = arg - 1;
	mpf_class arg_plus = arg + 1;
	mpf_class arg_minus_pow, arg_plus_pow;

	mpz_class two_in_pow, binom_prod;
	mpz_ui_pow_ui(two_in_pow.get_mpz_t(), 2, term);
	mpf_class coeff, step_res;

	for (std::size_t m = 0; m <= term; m++) {
		mpz_bin_uiui(binom1.get_mpz_t(), term + alpha, m);
		mpz_bin_uiui(binom2.get_mpz_t(), term + beta, term - m);
		mpz_mul(binom_prod.get_mpz_t(), binom1.get_mpz_t(), binom2.get_mpz_t());
		mpf_div(coeff.get_mpf_t(), ((mpf_class)binom_prod).get_mpf_t(), ((mpf_class)two_in_pow).get_mpf_t());
		mpf_pow_ui(arg_minus_pow.get_mpf_t(), arg_minus.get_mpf_t(), term - m);
		mpf_pow_ui(arg_plus_pow.get_mpf_t(), arg_plus.get_mpf_t(), m);
		step_res = mpf_class(0);
		mpf_mul(step_res.get_mpf_t(), coeff.get_mpf_t(), arg_minus_pow.get_mpf_t());
		mpf_mul(step_res.get_mpf_t(), step_res.get_mpf_t(), arg_plus_pow.get_mpf_t());
		mpf_add(res.get_mpf_t(), res.get_mpf_t(), step_res.get_mpf_t());
	}
	return res;
}

double Math::yacobi_polinom_init (std::size_t term, double arg, std::size_t alpha, std::size_t beta)
{
	double res = 0;
	for (std::size_t m = 0; m <= term; m++) {
		std::size_t binom1 = Math::binomial(term + alpha, m);
		std::size_t binom2 = Math::binomial(term + beta, term - m);
		std::size_t binom_prod = Math::product(binom1,binom2);
		std::size_t pow_of_two = Math::two_in_pow(term);
		double tmp_res = Math::divide(binom_prod, pow_of_two);
		res += tmp_res * std::pow(arg - 1, term - m) * std::pow(arg + 1, m);
	}
	return res;
}

double Math::yacobi_polinom_recur (std::size_t term, double arg, std::size_t alpha, std::size_t beta)
{
	if (alpha != 1 || beta != 0) {
		std::string text = "Recurent version doesn't support this parameters";
		throw std::invalid_argument(text);
	}

	if (term == 0) return 0;

	double res = 0;
	std::size_t two_pow = Math::two_in_pow(term);

	for (std::size_t m = 0; m <= term; m++) {
		std::size_t bprod = Math::binom_prod(term, m);
		double coeff = Math::divide(bprod, two_pow);
		res += coeff * std::pow(arg - 1, term - m) * std::pow(arg + 1, m);
	}

	return res;
	
}

std::size_t Math::next_prime(std::size_t prime, std::size_t max)
{
	if (prime == max) return prime+1;
	for (std::size_t next = prime+1; next <= max; next++) {
		bool is_prime = true;
		for (std::size_t denumer = 2; denumer < next; denumer++) {
			bool is_divide = (next % denumer == 0);
			bool is_equal = (next == denumer);
			if (is_divide && !is_equal) {
				is_prime = false;
				break;
			}
		}
		if (is_prime) return next;
	} 

	throw std::logic_error("I am your father!");
}

void Math::divide3to1 (std::size_t &x, std::size_t &y, std::size_t &z, std::size_t k)
{
	std::size_t prime = 2;
	while (k != 1) {
		if (k % prime == 0) {
			k /= prime;
			if (x % prime == 0) x /= prime;
			else if (y % prime == 0) y /= prime; 
			else if (z % prime == 0) z /= prime;
		} else prime = Math::next_prime(prime, k);
	}
}

std::size_t Math::binom_prod (std::size_t n, std::size_t m)
{
	if (n == m ) return n + 1;

	std::size_t res = 1;

	for (std::size_t k = 1; k <= m; k++) {

		/* res *= (n - k + 1) * (n - k + 2);
		res /= k * k; */

		std::size_t first = n - k + 1;
		std::size_t second = n - k + 2;

		if ( res % k == 0 ) res /= k;
		else if ( second % k == 0 ) second /= k;
		else if ( first % k == 0 ) first /= k;
		else Math::divide3to1(res, first, second, k);

		if ( res % k == 0 ) res /= k;
		else if ( second % k == 0 ) second /= k;
		else if ( first % k == 0 ) first /= k;
		else Math::divide3to1(res, first, second, k);

		res *= first * second;
	}

	return res;
}

float Math::inv_sqrt (float x)
{
	float xhalf = 0.5f * x;
	int i = *(int*)&x;
	i = 0x5f3759df - (i >> 1); 
	x = *(float*)&i;
	x = x*(1.5f-(xhalf*x*x));
	return x;
}

double Math::derivat3 ( std::function<double(double)> f, double x )
{
	double h = DERIVATIVE_STEP;
	return (f(x+h) - f(x-h)) / (2*h);
}

double Math::derivat4 ( std::function<double(double)> f, double x )
{
	double h = DERIVATIVE_STEP;
	return ( f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h) ) / (12*h);
}

double Math::lommel (std::size_t terms, int order, double W, double Z)
{
	if (Z < 1e-10) throw std::invalid_argument("Z=0 is not allowed in Math::lommel(term,n,W,Z)");

	const auto minus_one = mpf_class(-1).get_mpf_t();
	mpf_class res = mpf_class(0);
	mpf_class step_coeff, step_bessel, step_res;
	mpf_class gmp_W = W, gmp_Z = Z;

	for (std::size_t m = 0; m < terms; m++) {
		mpf_pow_ui(step_res.get_mpf_t(), minus_one, m);
		
		mpf_div(step_coeff.get_mpf_t(), gmp_W.get_mpf_t(), gmp_Z.get_mpf_t());
		mpf_pow_ui(step_coeff.get_mpf_t(), step_coeff.get_mpf_t(), order+2*m);
		mpf_mul(step_res.get_mpf_t(),step_res.get_mpf_t(),step_coeff.get_mpf_t());
		
		step_bessel = mpf_class(jn(order+2*m, Z));
		mpf_mul(step_res.get_mpf_t(),step_res.get_mpf_t(),step_bessel.get_mpf_t());
		
		mpf_add(res.get_mpf_t(),res.get_mpf_t(),step_res.get_mpf_t());
	}

	return res.get_d();
}
