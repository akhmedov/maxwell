#include "test.hpp"

using namespace std;

bool Test::next_prime()
{
	const vector<size_t> prime_list = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 
		41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
		127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 
		199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
		283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 
		383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
		467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 
		577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 
		661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 
		769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
		877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 
		983, 991, 997 };

	for (size_t iter = 0; iter < prime_list.size() - 2; iter++ ) {
		size_t next_prime = Math::next_prime(prime_list[iter], prime_list[iter+1]);
		if ( prime_list[iter + 1] != next_prime)
			return false;
	}
		
	for (size_t iter = 0; iter < prime_list.size() - 2; iter++ ) {
		if ( prime_list[iter + 1] != Math::next_prime(prime_list[iter], 997) )
			return false;
	}

	for (size_t iter = 0; iter < prime_list.size() - 2; iter++ ) {
		size_t next_prime = Math::next_prime(prime_list[iter], prime_list[iter+1]+1);
		if ( prime_list[iter + 1] != next_prime)
			return false;
	}

	/* // find max prime number for X min
	size_t max_num = std::numeric_limits<size_t>::max();
	size_t max_prime = 1;
	while (max_prime < max_num - 100)
		max_prime = Math::next_prime(max_prime, max_num);
	cout << max_prime << endl; */

	return true;
}

bool Test::binom_prod ()
{
	for (size_t n = 1; n <= 30; n++) {
		for (size_t m = 0; m <= n; m++) {
			size_t C1 = Math::binomial(n+1, m);
			size_t C2 = Math::binomial(n, n-m);
			size_t C12 = Math::binom_prod(n, m);
			if (C12 != C1*C2) return false;
		}
	}

	return true;
}

bool Test::yacobi_pol_property ()
{
	double eps = 10e-6;

	double table_value = Math::binomial(5,4);
	double value_basic = Math::yacobi_polinom_init(4, 1);
	double value_fast = Math::yacobi_polinom(4, 1).get_d();

	if (abs(table_value - value_basic) > eps) return false;
	if (abs(table_value - value_fast) > eps) return false;

	table_value = pow(-1,4) * Math::binomial(4, 4);
	value_basic = Math::yacobi_polinom_init(4, -1);
	value_fast = Math::yacobi_polinom(4, -1).get_d();

	if (abs(table_value - value_basic) > eps) return false;
	if (abs(table_value - value_fast) > eps) return false;

	size_t order = 4;
	size_t arg = 5;
	double res_basic = Math::yacobi_polinom_init(order, arg);
	double res_fast = Math::yacobi_polinom(order, arg).get_d();

	if (std::abs(res_basic - res_fast) > eps) return false;

	for (double arg = -1; arg < 1; arg += 0.01) {
		double res_fast = Math::yacobi_polinom(0, arg).get_d();
		if (std::abs(res_fast - 1) > eps) return false;
	}

	return true;
}

bool Test::yacobi_pol_compare ()
{
	double eps = 10e-4;

	for (size_t order = 1; order < 18; order++) {
		for (double arg = -3; arg < 3; arg += 0.01) {
			double res_basic = Math::yacobi_polinom_init(order, arg);
			double res_fast = Math::yacobi_polinom(order, arg).get_d();

			if (std::abs(res_basic - res_fast) > eps) return false;
		}		
	}

	return true;
}

bool Test::bessel_perp ()
{
	auto f = [] (double x) { return j0(x); } ;
	auto f_perp = [] (double x) { return -j1(x); } ;

	size_t iterator = 0;
	double error = 0;
	for (double arg = -5; arg < 5; arg += 10e-5) {
		error = abs( f_perp(arg) - Math::derivative(f,arg) );
		iterator++;
	}

	if (error/iterator > 10e-3) return false;
	else return true;
}

bool Test::simpson_I1 ()
{

	for (double a = 0.1; a <= 0.6; a += 0.5) {
		for (double b = 0.1; b <= 0.6; b += 0.5) {
			for (double c = 0.1; c <= 0.6; c += 0.5) {

				auto f = [a, b, c] (double x) {
					if (b == 0) return a * 0.5 * j1(a*x) * j0(c*x);
					if (x == 0) return 0.0;
					return a * j1(a*x) * j1(b*x) * j0(c*x) / (x * b);
				};

				double anal = MissileField::int_bessel_011(c, b, a);
				size_t bais = 10e5;
				Simpson I = Simpson(10*bais);
				double numerics = I.value(0, bais, f);

				double diff = abs(numerics - anal);
				double avarage = (numerics + anal) / 2;
				double error = 100 * diff / (2 * avarage);

				if (error > 1 && anal > 10e-10) return false;
			}
		}
	}

	return true;
}

bool Test::simpson_I2 ()
{
	for (double a = 0.1; a <= 0.6; a += 0.5) {
		for (double b = 0.1; b <= 0.6; b += 0.5) {
			for (double c = 0.1; c <= 0.6; c += 0.5) {

				auto f = [a, b, c] (double x) { 
					return c * j0(a*x) * j0(b*x) * j1(c*x);
				};

				double anal = MissileField::int_bessel_001(a, b, c);
				size_t bais = 10e5;
				Simpson I = Simpson(10*bais);
				double numerics = I.value(0, bais, f);

				double diff = abs(numerics - anal);
				double avarage = (numerics + anal) / 2;
				double error = 100 * diff / (2 * avarage);

				if (error > 1 && anal > 10e-10) return false;
			}
		}
	}

	return true;
}

bool Test::invers_sqrt ()
{
	double val = 2;
	double sqr = Math::inv_sqrt (val);
	if ( abs(sqr - 0.707) > 10e-5 ) return false;
	return true;
}

bool Test::field_getters ()
{
	double radius = 3.14159;
	double magnitude = 2.7;
	double eps_r = 2.1;
	double mu_r = 2.2;
	double chi = 10e-6;
	double sigma = 10e-7;

	UniformPlainDisk* source = new UniformPlainDisk(radius, magnitude);
	KerrMedium* medium = new KerrMedium(mu_r, eps_r, chi, sigma);
	MissileField* linear = new MissileField(source, medium);
	Test* non_linear = (Test*) new KerrAmendment(linear, medium);

	double permittivity = non_linear->nl_medium->relative_permittivity(0,0);
	double permeability = non_linear->nl_medium->relative_permeability(0,0);
	
	if ( abs(permittivity - eps_r) > 10e-16 ) return false;
	if ( abs(permeability - mu_r) > 10e-16 ) return false;

	if ( abs(non_linear->A0 - magnitude) > 10e-16 ) return false;
	if ( abs(non_linear->R - radius) > 10e-16 ) return false;

	return true;
}

bool Test::real_convolution ()
{
	double radius = 1;
	double magnitude = 2.7;
	double eps_r = 1;
	double mu_r = 1;
	double chi = 10e-6;
	double sigma = 10e-7;

	UniformPlainDisk* source = new UniformPlainDisk(radius, magnitude);
	KerrMedium* medium = new KerrMedium(mu_r, eps_r, chi, sigma);
	MissileField* linear = new MissileField(source, medium);
	Test* non_linear = (Test*) new KerrAmendment(linear, medium);

	complex<double> Erho = non_linear->electric_rho (0.2, 0, 0, 0.1);
	cout << "Result: " << Erho << endl;
	
	return false;
}

int main()
{
	cout << boolalpha;

	/* cout << "Math::Test::next_prime \t\t\t"; 
	cout.flush();
	cout << Test::next_prime() << endl;

	cout << "Math::Test::binom_prod \t\t\t"; 
	cout.flush();
	cout << Test::binom_prod() << endl;

	cout << "Math::Test::yacobi_pol_property \t"; 
	cout.flush();
	cout << Test::yacobi_pol_property() << endl;

	cout << "Math::Test::yacobi_pol_compare \t\t"; 
	cout.flush();
	cout << Test::yacobi_pol_compare() << endl;

	cout << "Math::Test::invers_sqrt \t\t"; 
	cout.flush();
	cout << Test::invers_sqrt() << endl;

	cout << "Math::Test::devidative \t\t\t"; 
	cout.flush();
	cout << Test::bessel_perp() << endl;

	cout << "MissileField::Test::simpson_I1 \t\t"; 
	cout.flush();
	cout << Test::simpson_I1() << endl;

	cout << "MissileField::Test::simpson_I2 \t\t"; 
	cout.flush();
	cout << Test::simpson_I2() << endl;

	cout << "KerrAmendment::Test::field_getters \t"; 
	cout.flush();
	cout << Test::field_getters() << endl; */

	cout << "KerrAmendment::Test::real_convolution \t";
	cout.flush();
	cout << Test::real_convolution() << endl;
	
	return 0;
}
