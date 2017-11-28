//
//  integral.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 17.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"

// =========================================================================

Simpson::Simpson( std::size_t terms )
: quadr_terms(terms) { }

double Simpson::value (double from, double to, std::function<double(double)> func) const
{
	if (from >= to) throw std::invalid_argument("Low bound of integral is bigger then upper.");

	double h = (to - from) / this->quadr_terms;
	double I = 0;

	/* for (std::size_t i = 0; i < terms; i++) {
		double a = from + h * i;
		double b = from + h * (i+1);
		I += func(a) + 4 * func((a + b)/2) + func(b);
	} */

	double a = from;
	double b = from + h;
	double f_a = func(a);
	double f_ab = func((a + b)/2);
	double f_b = func(b);

	while (b <= to) {
		I += f_a + 4 * f_ab + f_b;
		a = b;
		b += h;
		f_a = f_b;
		f_ab = func((a + b)/2);
		f_b = func(b);
	}

	return I * (to - from) / this->quadr_terms / 6;
}

// =========================================================================

MonteCarlo::MonteCarlo (std::size_t rolls, const std::vector<std::pair<double,double>> &limits )
{ 
	this->rand_rolls = rolls;
	this->volume = 1;
	
	for (auto l : limits) {
		this->volume *= std::abs(l.second - l.first);

		std::uniform_real_distribution<double> uniform(l.first, l.second);
		this->distribution.push_back(uniform);

		std::random_device rdev;
		this->generator.push_back(std::mt19937_64(rdev()));
	}
}

std::valarray<double> MonteCarlo::random_array ()
{
	// TODO: experement with perfomrnce - try pair container and auto iterator
	std::valarray<double> rand(this->generator.size());
	for (std::size_t i = 0; i < this->generator.size(); i++) {
		double value = this->distribution[i](this->generator[i]);
		rand[i] = value;
	}
	return rand;
}

double MonteCarlo::value ( const std::function<double(double)> &func )
{
	if (this->distribution.size() != 1) throw std::invalid_argument("Wrong dimention for Monte-Carlo!");
	
	double res = 0;
	std::size_t current_roll = 0;
	while (current_roll <= rand_rolls) {
		std::valarray<double> rand = random_array();
		res += func(rand[0]);
		current_roll++;
	}

	return this->volume * res / rand_rolls;	
}

double MonteCarlo::value ( const std::function<double(double, double)> &func )
{
	if (this->distribution.size() != 2) throw std::invalid_argument("Wrong dimention for Monte-Carlo!");

	double res = 0;
	std::size_t current_roll = 0;
	while (current_roll <= rand_rolls) {
		std::valarray<double> rand = random_array();
		res += func(rand[0], rand[1]);
		current_roll++;
	}

	return this->volume * res / rand_rolls;
}

double MonteCarlo::value ( const std::function<double(double, double, double)> &func )
{
	if (this->distribution.size() != 3) throw std::invalid_argument("Wrong dimention for Monte-Carlo!");

	double res = 0;
	std::size_t current_roll = 0;
	while (current_roll <= rand_rolls) {
		std::valarray<double> rand = random_array();
		res += func(rand[0], rand[1], rand[2]);
		current_roll++;
	}

	return this->volume * res / rand_rolls;
}

double MonteCarlo::value ( const std::function<double(double, double, double, double)> &func )
{
	if (this->distribution.size() != 4) throw std::invalid_argument("Wrong dimention for Monte-Carlo!");
	
	double res = 0;
	std::size_t current_roll = 0;
	while (current_roll <= rand_rolls) {
		// if () std::cout <<  << std::endl;
		std::valarray<double> rand = random_array();
		res += func(rand[0], rand[1], rand[2], rand[3]);
		current_roll++;
	}

	return this->volume * res / rand_rolls;
}

// =========================================================================

GaussLaguerre::GaussLaguerre ()
: quadr_terms(polynom_data.size()) { }

GaussLaguerre::GaussLaguerre (std::size_t terms)
: quadr_terms(polynom_data.size()) {
	if (terms > this->polynom_data.size()) 
		throw std::logic_error("Max terms number extended!");
}

double GaussLaguerre::value (std::function<double(double)> func) const
{
	double sum = 0;
	for (auto i : this->polynom_data)
		sum += i.second * func(i.first) * std::exp(i.first);
	return sum;
}

// =========================================================================

/*

const static short unsigned numberOfStepsInPeriod = 7;

double BesselIntegral::simpson(const WorldPoint& event) const
{
	double maxCoeff = 1, minCoeff = besselCoefficient[0].second;
	for (auto i : besselCoefficient) {
		if (maxCoeff < i.second) maxCoeff = i.second;
		if (minCoeff > i.second && i.second != 0) minCoeff = i.second;
	}
	
	double simpsonStep = 2 * M_PI / numberOfStepsInPeriod / maxCoeff;
	double properIntervel = (1 + 0.5) * 2 * M_PI / minCoeff;
	unsigned stepCount = properIntervel / simpsonStep;
	
	double res = 0, prevRes = 2 * accurenceForSimpson;
	unsigned iter = 0;
	while (std::abs(res - prevRes) > accurenceForSimpson) {
		prevRes = res;
		for (unsigned i = 0; i <= stepCount; i++) {
			double a = iter * properIntervel + i * simpsonStep;
			double b = iter * properIntervel + (i + 1) * simpsonStep;
			double intValue = underIntegral(a, event);
			intValue += 4 * underIntegral((a+b)/2, event);
			intValue += underIntegral(b, event);
			intValue *= b - a;
			// if (std::isnan(intValue)) 
			//         "Warning! Nan!" << std::endl;
			// else 
			res += intValue;
		}
		iter++;
	}
	
	return res / 6;
}

*/
