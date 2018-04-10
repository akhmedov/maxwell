//
//  unit_test.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 27.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "unit_test.hpp"

using namespace std;

bool UnitTest::next_prime()
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

bool UnitTest::binom_prod ()
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

bool UnitTest::yacobi_pol_property ()
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

bool UnitTest::yacobi_pol_compare ()
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

bool UnitTest::bessel_perp ()
{
	auto f = [] (double x) { return j0(x); } ;
	auto f_perp = [] (double x) { return -j1(x); } ;

	size_t iterator = 0;
	double error = 0;
	for (double arg = -5; arg < 5; arg += 10e-5) {
		error = abs( f_perp(arg) - Math::derivat4(f,arg) );
		iterator++;
	}

	if (error/iterator > 10e-3) return false;
	else return true;
}

bool UnitTest::simpson_I1 ()
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

bool UnitTest::simpson_I2 ()
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

bool UnitTest::invers_sqrt ()
{
	double val = 2;
	double sqr = Math::inv_sqrt (val);
	if ( abs(sqr - 0.707) > 10e-5 ) return false;
	return true;
}

bool UnitTest::field_getters ()
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
	UnitTest* non_linear = (UnitTest*) new KerrAmendment(linear, medium, source);

	double permittivity = non_linear->nl_medium->relative_permittivity(0,0);
	double permeability = non_linear->nl_medium->relative_permeability(0,0);
	
	if ( abs(permittivity - eps_r) > 10e-16 ) return false;
	if ( abs(permeability - mu_r) > 10e-16 ) return false;

	if ( abs(non_linear->A0 - magnitude) > 10e-16 ) return false;
	if ( abs(non_linear->R - radius) > 10e-16 ) return false;

	return true;
}

bool UnitTest::monte_carlo_vector ()
{
	vector<pair<double,double>> limits;
	limits.push_back(make_pair(0,1));
	limits.push_back(make_pair(0,1));

	MonteCarlo distr = MonteCarlo (3000, limits );
	std::valarray<double> val = distr.random_array();

	return false;
}

bool UnitTest::monte_carlo_integral ()
{
	double radius = 1;

	auto func = [] (double rho, double phi, double theta) {
		UNUSED(phi);
		return rho * rho * sin(theta);
	};

	auto int_func = [radius] () {
		double R3 = radius * radius * radius;
		return 4 * M_PI * R3 / 3;
	};

	vector<pair<double,double>> limits;
	limits.push_back( make_pair(0,radius) );
	limits.push_back( make_pair(0,M_PI_2) );
	limits.push_back( make_pair(0,M_PI_2) );

	MonteCarlo integral = MonteCarlo (10e8, limits );
	double volume = 8 * integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 1 ? true : false;
}

bool UnitTest::monte_carlo_improper ()
{
	double a = 1, b = 1;

	auto func = [a,b] (double x, double y, double z) {
		return exp(-y-z) * j0(a*sqrt(x*y)) * j0(b*sqrt(x*z));
	};

	auto int_func = [a,b] () {
		return a*a/4 + b*b/4;
	};

	vector<pair<double,double>> limits;
	limits.push_back(make_pair(0,10e3));
	limits.push_back(make_pair(0,10e3));
	limits.push_back(make_pair(0,10e3));

	MonteCarlo integral = MonteCarlo(10e8, limits);
	double volume = integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 10 ? true : false;
}

bool UnitTest::imptoper_int_bessel ()
{
	double a = 2;

	auto func = [a] (double x) {
		return j0(a * x);
	};

	auto int_func = [a] () {
		return 1 / a;
	};

	vector<pair<double,double>> limits;
	limits.push_back(make_pair(0,10e3));
	MonteCarlo integral = MonteCarlo(10e8, limits);
	double monte_carlo = integral.value(func);
	double error_mc = 100 * abs(monte_carlo-int_func()) / int_func();

	size_t bais = 10e5;
	Simpson I = Simpson(10*bais);
	double simpson = I.value(0, bais, func);
	double error_s = 100 * abs(simpson-int_func()) / int_func();

	return (error_s < 1 ? true : false) && (error_mc < 4 ? true : false);
}

bool UnitTest::simpson2d ()
{
	double radius = 1;

	auto func = [] (double rho, double phi) {
		UNUSED(phi);
		return rho;
	};

	double etalon = M_PI * radius * radius;

	vector<tuple<double,size_t,double>> limits;
	limits.push_back( make_tuple(0,300,radius) );
	limits.push_back( make_tuple(0,300,M_PI_2) );

	Simpson2D integral = Simpson2D(limits);
	double square = 4 * integral.value(func);
	double error = 100 * abs(square-etalon) / etalon;

	return error < 1 ? true : false;
}

bool UnitTest::simpson2d_line ()
{
	double a = 2;
	auto func = [] (double x, double y) { return x * y; };

	Simpson2D_line integral = Simpson2D_line();
	integral.first_limit(0, 400, a);
	auto min_y  = [] (double x) { UNUSED(x); return 0;     };
	auto term_y = [] (double x) { return (std::size_t) 200*x; };
	auto max_y  = [] (double x) { return x;     };
	integral.second_limit(min_y, term_y, max_y);

	double res = integral.value(func);
	double error = 100 * abs(res-2) / 2;
	return error < 1 ? true : false;
}

bool UnitTest::simpson2d_line2 ()
{
	double radius = 1;

	auto func = [] (double rho, double phi) {
		UNUSED(phi);
		return rho;
	};

	Simpson2D_line integral = Simpson2D_line();
	integral.first_limit(0, 300, radius);
	auto min_y  = [] (double x) { UNUSED(x); return 0;      };
	auto term_y = [] (double x) { UNUSED(x); return 300;    };
	auto max_y  = [] (double x) { UNUSED(x); return M_PI_2; };
	integral.second_limit(min_y, term_y, max_y);

	double res = 4 * integral.value(func);
	double error = 100 * abs(res-M_PI) / M_PI;
	return error < 1 ? true : false;
}

bool UnitTest::simpson3d ()
{
	double radius = 1;

	auto func = [] (double rho, double phi, double theta) {
		UNUSED(phi);
		return rho * rho * sin(theta);
	};

	auto int_func = [radius] () {
		double R3 = radius * radius * radius;
		return 4 * M_PI * R3 / 3;
	};

	vector<tuple<double,size_t,double>> limits;
	limits.push_back( make_tuple(0,200,radius) );
	limits.push_back( make_tuple(0,200,M_PI_2) );
	limits.push_back( make_tuple(0,200,M_PI_2) );

	Simpson3D integral = Simpson3D(limits);
	double volume = 8 * integral.value(func);
	double error = 100 * abs(volume-int_func()) / int_func();

	return error < 1 ? true : false;
}

bool UnitTest::real_convolution ()
{
	double radius = 1;
	double magnitude = 100;
	double eps_r = 1;
	double mu_r = 1;
	double chi = 10e-15;
	double sigma = 0;

	UniformPlainDisk* source = new UniformPlainDisk(radius, magnitude);
	KerrMedium* medium = new KerrMedium(mu_r, eps_r, chi, sigma);
	MissileField* linear = new MissileField(source, medium);
	UnitTest* non_linear = (UnitTest*) new KerrAmendment(linear, medium, source);

	double nonlinear_Erho = non_linear->electric_rho (2.1, 0.8, 0, 2);
	double linear_Erho = linear->electric_rho (2.1, 0.8, 0, 2);
	cout << "Result: " << linear_Erho << ' ' << nonlinear_Erho << endl;
	
	return false;
}

bool UnitTest::I1_time_partder ()
{
	for (double rho = 0; rho <= 1.3; rho += 0.25) {
		for (double R = 0.1; R <= 1.3; R += 0.2) {
			for (double ct = 0.1; ct <= 1.4; ct += 0.3) {
				for (double z = 1e-2; z <= ct - 0.05; z += 0.3) {

					auto I1 = [rho, R, z] (double ct) {
						double vt_z = std::sqrt(ct * ct - z * z);
						return MissileField::int_bessel_011(vt_z,rho,R);
					};
					
					double anal = KerrAmendment::int_bessel_011_perp(ct,z,rho,R);
					auto num3 = Math::derivat3(I1, ct);
					auto num4 = Math::derivat4(I1, ct);

					double error3 = anal ? std::abs(anal - num3) / anal : 0;
					double error4 = anal ? std::abs(anal - num4) / anal : 0;
					if (100*error3 > 10) return false;
					if (100*error4 > 5) return false;
				}
			}
		}
	}

	return true;
}

bool UnitTest::I2_time_partder ()
{
	for (double rho = 0; rho <= 1.3; rho += 0.25) {
		for (double R = 0.1; R <= 1.3; R += 0.2) {
			for (double ct = 0.1; ct <= 1.4; ct += 0.3) {
				for (double z = 1e-2; z <= ct - 0.05; z += 0.3) {

					auto I2 = [rho, R, z] (double ct) {
						double vt_z = std::sqrt(ct * ct - z * z);
						return MissileField::int_bessel_001(vt_z,rho,R); 
					};
					
					double anal = KerrAmendment::int_bessel_001_perp(ct,z,rho,R);
					auto num3 = Math::derivat3(I2, ct);
					auto num4 = Math::derivat4(I2, ct);

					double error3 = anal ? std::abs(anal - num3) / anal : 0;
					double error4 = anal ? std::abs(anal - num4) / anal : 0;
					if (100*error3 > 10) return false;
					if (100*error4 > 5) return false;
				}
			}
		}
	}

	return true;
}

bool UnitTest::simpson_runge ()
{
	double accuracy = 1; // %

	double from = 0;
	double to = M_PI_2;

	auto f = [] (double x) {
		return std::sin(x);
	};

	double I = SimpsonRunge(1, 1).value(from, to, f);
	return 100*std::abs(I-1) < accuracy ? true : false;
}

bool UnitTest::simpson_runge_2d ()
{
	size_t init_units = 1;
	double error = 2; /* % */

	auto f = [] (double x, double y) {
		return 4 * x * y;
	};

	SimpsonRunge integral = SimpsonRunge(init_units, error);

	double I = integral.value(0, 1,
		[f, &integral] (double x) {
			return integral.value(0, 1, [f,x] (double y) {
				return f(x,y); 
			} );
		} );

	return 100 * abs(I - 1) < 2*error ? true : false;
}

int main()
{
	cout << boolalpha;

	cout << "Math::UnitTest::next_prime \t\t\t"; 
	cout.flush();
	cout << UnitTest::next_prime() << endl;

	cout << "Math::UnitTest::binom_prod \t\t\t"; 
	cout.flush();
	cout << UnitTest::binom_prod() << endl;

	cout << "Math::UnitTest::yacobi_pol_property \t\t"; 
	cout.flush();
	cout << UnitTest::yacobi_pol_property() << endl;

	cout << "Math::UnitTest::yacobi_pol_compare \t\t"; 
	cout.flush();
	cout << UnitTest::yacobi_pol_compare() << endl;

	/* cout << "Math::UnitTest::invers_sqrt \t\t\t"; 
	cout.flush();
	cout << UnitTest::invers_sqrt() << endl; */

	cout << "MissileField::UnitTest::simpson_runge \t\t"; 
	cout.flush();
	cout << UnitTest::simpson_runge() << endl;

	cout << "MissileField::UnitTest::simpson_runge_2d \t"; 
	cout.flush();
	cout << UnitTest::simpson_runge_2d() << endl;

	cout << "MissileField::UnitTest::simpson_I1 \t\t"; 
	cout.flush();
	cout << UnitTest::simpson_I1() << endl;

	cout << "MissileField::UnitTest::simpson_I2 \t\t"; 
	cout.flush();
	cout << UnitTest::simpson_I2() << endl;

	cout << "Math::UnitTest::simpson2d \t\t\t"; 
	cout.flush();
	cout << UnitTest::simpson2d() << endl;

	cout << "Math::UnitTest::simpson2d_line \t\t\t"; 
	cout.flush();
	cout << UnitTest::simpson2d_line() << endl;

	cout << "Math::UnitTest::simpson2d_line2 \t\t"; 
	cout.flush();
	cout << UnitTest::simpson2d_line2() << endl;

	cout << "Math::UnitTest::simpson3d \t\t\t"; 
	cout.flush();
	cout << UnitTest::simpson3d() << endl;

	cout << "KerrAmendment::UnitTest::field_getters \t\t"; 
	cout.flush();
	cout << UnitTest::field_getters() << endl;

	cout << "Math::UnitTest::deridative \t\t\t";
	cout.flush();
	cout << UnitTest::bessel_perp() << endl;

	cout << "KerrAmendment::UnitTest::I1_time_partder \t"; 
	cout.flush();
	cout << UnitTest::I1_time_partder() << endl;

	cout << "KerrAmendment::UnitTest::I2_time_partder \t"; 
	cout.flush();
	cout << UnitTest::I2_time_partder() << endl;

	/* cout << "Math::UnitTest::monte_carlo_integral \t";
	cout.flush();
	cout << UnitTest::monte_carlo_integral() << endl;

	cout << "Math::UnitTest::monte_carlo_vector \t\t"; 
	cout.flush();
	cout << UnitTest::monte_carlo_vector() << endl;

	cout << "Math::UnitTest::monte_carlo_improper \t";
	cout.flush();
	cout << UnitTest::monte_carlo_improper() << endl;

	cout << "Math::UnitTest::imptoper_int_bessel \t";
	cout.flush();
	cout << UnitTest::imptoper_int_bessel() << endl; */

	/* cout << "KerrAmendment::UnitTest::real_convolution \t";
	cout.flush();
	cout << UnitTest::real_convolution() << endl; */

	return 0;
}
