//
//  integral.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 17.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include <functional>
#include <random>
#include <valarray>
#include <vector>

using vector_pair_dd = std::vector<std::pair<double,double>>;
using vector_tuple_did = std::vector<std::tuple<double,std::size_t,double>> ;

struct Integral {}; // TODO: Why it's needed?

struct Simpson : public Integral {
	Simpson(std::size_t terms) : quadr_terms(terms) {}
	double value (double from, double to, const std::function<double(double)> &f) const; // TODO: Use aliases for all std::funciton types?

private:
	std::size_t quadr_terms{};
};

struct SimpsonRunge : public Integral {
	SimpsonRunge (std::size_t node_number, double precision /* % */, std::size_t nodes = 1e4)
		: init_nodes(node_number), max_nodes(nodes), epsilon(precision)
	{
		if (node_number == 0) throw std::invalid_argument("Zero node number is not allowed.");
	}
	
	double value (double from, double to, const std::function<double(double)> &f);
	std::size_t current_units () const { return running_units; }

private:
	std::size_t init_nodes{};
	std::size_t max_nodes{};
	double epsilon{};
	std::size_t running_units{ 1 };
};

struct Simpson2D : public Integral {
	Simpson2D ( const vector_tuple_did &limits )
		: x_min(std::get<0>(limits[0])), x_terms(std::get<1>(limits[0])), x_max(std::get<2>(limits[0])),
  		  y_min(std::get<0>(limits[1])), y_terms(std::get<1>(limits[1])), y_max(std::get<2>(limits[1])) 
	{ 
		if (limits.size() != 2) 
			throw std::invalid_argument("Simpson2D: Only 3dim is implemented!");
		if ( (x_min >= x_max) || (y_min >= y_max) )
			throw std::invalid_argument("Low bound of integral is bigger then upper.");
	}

	double value (const std::function<double(double,double)> &func) const;

private:
	double x_min{};
	std::size_t x_terms{};
	double x_max{};

	double y_min{};
	std::size_t y_terms{};
	double y_max{};
};

struct Simpson2D_line : public Integral {
	void first_limit (double from, std::size_t terms, double to)
	{
		x_min = from;
		x_terms = terms;
		x_max = to;
	}

	void second_limit (const std::function<double(double)> &from, const std::function<std::size_t(double)> &terms, const std::function<double(double)> &to)
	{
		y_min = from;
		y_terms = terms;
		y_max = to;
	}

	double value (const std::function<double(double,double)> &func) const;

private:
	double x_min{};
	std::size_t x_terms{};
	double x_max{};

	std::function<double(double)> y_min;
	std::function<std::size_t(double)> y_terms;
	std::function<double(double)> y_max;
};

struct Simpson3D : public Integral {
	Simpson3D ( const vector_tuple_did &limits );
	double value (const std::function<double(double,double,double)> &func) const;

private:
	double x_min{};
	std::size_t x_terms{};
	double x_max{};

	double y_min{};
	std::size_t y_terms{};
	double y_max{};

	double z_min{};
	std::size_t z_terms{};
	double z_max{};
};

struct MonteCarlo : public Integral {
	MonteCarlo ( std::size_t rolls, const vector_pair_dd &limits );
	std::valarray<double> random_array ();

	double value ( const std::function<double(double)> &func );
	double value ( const std::function<double(double, double)> &func );
	double value ( const std::function<double(double, double, double)> &func );
	double value ( const std::function<double(double, double, double, double)> &func );

private:
	double volume{};
	std::size_t rand_rolls{};
	// TODO: implement std::normal_distribution<>
	std::vector<std::mt19937_64> generator;
	std::vector<std::uniform_real_distribution<double>> distribution;
};

struct GaussLaguerre : public Integral {
	GaussLaguerre () = default;
	GaussLaguerre (std::size_t terms)
	{
		if (terms > polynom_data.size()) 
			throw std::logic_error("Max terms number extended!");
	}

	double value (std::function<double(double)> f) const;
	
private:
	std::vector<std::pair<double, double>> polynom_data {
		/*				root				weight				*/
		std::make_pair( 0.0933078120170,	0.218234885940		),
		std::make_pair( 0.4926917403020,	0.342210177923		),
		std::make_pair( 1.2155954120710,	0.263027577942		),
		std::make_pair( 2.2699495262040,	0.126425818106		),
		std::make_pair( 3.6676227217510,	0.040206864921		),
		std::make_pair( 5.4253366274140,	0.856387780361e-2	),
		std::make_pair( 7.5659162266130,	0.121243614721e-2	),
		std::make_pair( 10.120228568019,	0.111674392344e-3	),
		std::make_pair( 13.130282482176,	0.645992676202e-5	),
		std::make_pair( 16.654407708330,	0.222631690710e-6	),
		std::make_pair( 20.776478899449,	0.422743038498e-8	),
		std::make_pair( 25.623894226729,	0.392189726704e-10	),
		std::make_pair( 31.407519169754,	0.145651526407e-12	),
		std::make_pair( 38.530683306486,	0.148302705111e-15	),
		std::make_pair( 48.026085572686,	0.160059490621e-19	)
	};
	std::size_t quadr_terms{ polynom_data.size() };
};
