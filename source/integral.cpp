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

SimpsonMultiDim::SimpsonMultiDim ( const vector_tuple_did &limits )
: w_min(std::get<0>(limits[0])), w_terms(std::get<1>(limits[0])), w_max(std::get<2>(limits[0])),
  x_min(std::get<0>(limits[1])), x_terms(std::get<1>(limits[1])), x_max(std::get<2>(limits[1])),
  y_min(std::get<0>(limits[2])), y_terms(std::get<1>(limits[2])), y_max(std::get<2>(limits[2])),
  z_min(std::get<0>(limits[3])), z_terms(std::get<1>(limits[3])), z_max(std::get<2>(limits[3])) 
{ 
	if (limits.size() != 4) 
		throw std::invalid_argument("SimpsonMultiDim: Only 4dim is implemented!");
	if ( (w_min >= w_max) || (x_min >= x_max) || (y_min >= y_max) || (z_min >= z_max) )
		throw std::invalid_argument("Low bound of integral is bigger then upper.");
}

double SimpsonMultiDim::value (const std::function<double(double,double,double,double)> &func) const
{
	double hw = (this->w_max - this->w_min) / this->w_terms;
	double hx = (this->x_max - this->x_min) / this->x_terms;
	double hy = (this->y_max - this->y_min) / this->y_terms;
	double hz = (this->z_max - this->z_min) / this->z_terms;

	double total_volume = (this->w_max-this->w_min) * (this->x_max-this->x_min) * (this->y_max-this->y_min) * (this->z_max-this->z_min);
	double total_terms = this->w_terms * this->x_terms * this->y_terms * this->z_terms;
	double I = 0;

	std::vector<std::vector<double>> 
	arg_grid(4, std::vector<double>
			(3, 0.0));

	arg_grid[0][0] = this->w_min; 							// w_a
	arg_grid[0][2] = this->w_min + hw;						// w_ab
	arg_grid[0][1] = (arg_grid[0][0] + arg_grid[0][2]) / 2; // w_b

	arg_grid[1][0] = this->x_min;							// x_a
	arg_grid[1][2] = this->x_min + hx;						// x_ab
	arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2;	// x_b

	arg_grid[2][0] = this->y_min;							// y_a
	arg_grid[2][2] = this->y_min + hy;						// y_ab
	arg_grid[2][1] = (arg_grid[2][0] + arg_grid[2][2]) / 2; // y b

	arg_grid[3][0] = this->z_min;							// z_a
	arg_grid[3][2] = this->z_min + hz;						// z_ab
	arg_grid[3][1] = (arg_grid[3][0] + arg_grid[3][2]) / 2; // z_b

	// init setup for fun_grid
	std::vector<std::vector<std::vector<std::vector<double>>>> 
	fun_grid(3, std::vector<std::vector<std::vector<double>>>
			(3, std::vector<std::vector<double>>
			(3, std::vector<double>
			(3, 0.0))));

	for (std::size_t w = 0; w < 3; w++) {
		for (std::size_t x = 0; x < 3; x++) {
			for (std::size_t y = 0; y < 3; y++) {
				for (std::size_t z = 0; z < 3; z++) {
					fun_grid[w][x][y][z] = func(arg_grid[0][w],
												arg_grid[1][x],
												arg_grid[2][y],
												arg_grid[3][z]);
				}
			}
		}
	}

	/* init setup for c_grid */
	std::vector<std::vector<std::vector<std::vector<double>>>> 
	c_grid  (3, std::vector<std::vector<std::vector<double>>>
			(3, std::vector<std::vector<double>>
			(3, std::vector<double>
			(3, 1.0))));

	for (std::size_t w = 0; w < 3; w++) {
		for (std::size_t x = 0; x < 3; x++) {
			for (std::size_t y = 0; y < 3; y++) {
				for (std::size_t z = 0; z < 3; z++) {
					if (w == 1) c_grid[w][x][y][z] *= 4.0;
					if (x == 1) c_grid[w][x][y][z] *= 4.0;
					if (y == 1) c_grid[w][x][y][z] *= 4.0;
					if (z == 1) c_grid[w][x][y][z] *= 4.0;
				}
			}
		}
	}

	/* 
	w_max    - 100
	arg_grid - x
	*/


	while (arg_grid[0][2] <= this->w_max) {
		std::cout << arg_grid[0][2] * this->w_max / 100 << '%' << std::endl;
		arg_grid[1][0] = this->x_min;
		arg_grid[1][2] = this->x_min + hx;
		arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2;
		
		while (arg_grid[1][2] <= this->x_max) {
			arg_grid[2][0] = this->y_min;
			arg_grid[2][2] = this->y_min + hy;
			arg_grid[2][1] = (arg_grid[2][0] + arg_grid[2][2]) / 2;
			
			while (arg_grid[2][2] <= this->y_max) {
				arg_grid[3][0] = this->z_min;
				arg_grid[3][2] = this->z_min + hz;
				arg_grid[3][1] = (arg_grid[3][0] + arg_grid[3][2]) / 2;
				
				while (arg_grid[3][2] <= this->z_max) {

					for (std::size_t w = 0; w < 3; w++) {
						for (std::size_t x = 0; x < 3; x++) {
							for (std::size_t y = 0; y < 3; y++) {
								for (std::size_t z = 0; z < 3; z++) {
									I += c_grid[w][x][y][z] * fun_grid[w][x][y][z];
								}
							}
						}
					}

					arg_grid[3][0] = arg_grid[3][2]; // UPD: z_min
					arg_grid[3][2] += hz; // UPD: z_max
					arg_grid[3][1] = (arg_grid[3][0] + arg_grid[3][2]) / 2; // UPD: z_mid
					
					// fuc_grid reset: TODO (not optimum)
					for (std::size_t w = 0; w < 3; w++) {
						for (std::size_t x = 0; x < 3; x++) {
							for (std::size_t y = 0; y < 3; y++) {
								for (std::size_t z = 0; z < 3; z++) {
									fun_grid[w][x][y][z] = func(arg_grid[0][w],
																arg_grid[1][x],
																arg_grid[2][y],
																arg_grid[3][z]);
								}
							}
						}
					}

				}
				arg_grid[2][0] = arg_grid[2][2]; // UPD: y_min
				arg_grid[2][2] += hy; // UPD: y_max
				arg_grid[2][1] = (arg_grid[2][0] + arg_grid[2][2]) / 2; // UPD: y_mid
			}
			arg_grid[1][0] = arg_grid[1][2]; // UPD: x_min
			arg_grid[1][2] += hx; // UPD: x_max
			arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2; // UPD: x_mid
		}
		arg_grid[0][0] = arg_grid[0][2]; // UPD: w_min
		arg_grid[0][2] += hw; // UPD: w_max
		arg_grid[0][1] = (arg_grid[0][0] + arg_grid[0][2]) / 2; // UPD: w_mid
	}

	return total_volume / total_terms / 1296.0 * I;
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
		std::valarray<double> rand = this->random_array();
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
		std::valarray<double> rand = this->random_array();
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
		std::valarray<double> rand = this->random_array();
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
		std::valarray<double> rand = this->random_array();
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
