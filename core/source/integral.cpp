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

double Simpson::value (double from, double to, const std::function<double(double)> &func) const
{
	if (from >= to) return 0;

	double h = std::abs(to - from) / this->quadr_terms;
	double I = 0;

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

SimpsonRunge::SimpsonRunge (std::size_t min_node_number, double accurence, std::size_t max_node_number)
: init_nodes(min_node_number), max_nodes(max_node_number), epsilon(accurence)
{ 
	if (min_node_number == 0) throw std::invalid_argument("Zero node number is not allowed.");
	this->running_units = 1;
}

double SimpsonRunge::value (double from, double to, const std::function<double(double)> &func)
{
	if (from >= to) return 0;
	if (std::isnan(from)) throw std::logic_error("FROM limit is nan.");
	if (std::isnan(to))   throw std::logic_error("TO limit is nan.");

	// First iteration with Simpson quadtature

	double newI = 0;
	std::pair<double,double> old_sum(0,0); // { f_a + f_b , f_ab }
	std::pair<double,double> new_sum(0,0); // { f_a + f_b , f_ab }
	double h = std::abs(to - from) / this->init_nodes;

	double a = from;
	double b = from + h;
	double f_a = func(a);
	double f_ab = func(a+h/2);
	double f_b = func(b);

	while (b <= to) {
		old_sum.first += f_a + f_b;
		old_sum.second += f_ab;
		a = b;
		b += h;
		f_a = f_b;
		f_ab = func(a+h/2);
		f_b = func(b);
	}

	double oldI = (old_sum.first + 4 * old_sum.second) / this->init_nodes;

	// Runge rule for Simpson quadtature for augumentation

	for (std::size_t nodes = 2 * this->init_nodes; nodes <= this->max_nodes; nodes *= 2) {

		this->running_units *= 2;

		h /= 2;
		a = from + h;
		new_sum.first = 0;
		new_sum.second = 0;

		while (a < to) {
			new_sum.first  += func(a) / 2; // TODO: why devided by 2 ???
			new_sum.second += func(a-h/2);
			a += h;
		}

		new_sum.first -= func(to);
		new_sum.first += old_sum.first + old_sum.second;
		newI = (new_sum.first + 4 * new_sum.second) / nodes; // TODO: (node + 1) ???

		// refactor of exit criteria
		if (std::abs(newI) < 1e-5 ) return 0;
		double error = 100 * std::abs((oldI - newI) / newI);
		if (error < this->epsilon)  {
			newI = 32 * newI / 31 - oldI / 31; // Runge error linearization
			return (to - from) / 6 * newI;
		}

		oldI = newI;
		old_sum.first = new_sum.first;
		old_sum.second = new_sum.second;
	}

	throw (to - from) / 6 * newI;
}

std::size_t SimpsonRunge::current_units () const
{
	return this->running_units;
}

// =========================================================================


Simpson2D_line::Simpson2D_line () { }

void Simpson2D_line::first_limit (double from, std::size_t terms, double to)
{
	this->x_min = from;
	this->x_terms = terms;
	this->x_max = to;
}

void Simpson2D_line::second_limit (const std::function<double(double)> &from, 
								   const std::function<std::size_t(double)> &terms, 
								   const std::function<double(double)> &to)
{
	this->y_min = from;
	this->y_terms = terms;
	this->y_max = to;
}

double Simpson2D_line::value (const std::function<double(double,double)> &func) const
{
	double I = 0;
	double x = this->x_min;
	double hx = (this->x_max - this->x_min) / this->x_terms;

	if (this->x_max <= this->x_min) return 0;

	while (x <= this->x_max - hx) {

		double y = this->y_min(x);
		double hy = (this->y_max(x) - this->y_min(x)) / this->y_terms(x);
		double Iy = 0;

		if (this->y_max(x) > this->y_min(x)) {

			while (y <= this->y_max(x) - hy) {
				Iy += func(x,y);			Iy += 4 * func(x+hx/2,y);		Iy += func(x+hx,y);
				Iy += 4 * func(x,y+hy/2);	Iy += 16 * func(x+hx/2,y+hy/2);	Iy += 4 * func(x+hx,y+hy/2);
				Iy += func(x,y+hy);			Iy += 4 * func(x+hx/2,y+hy);	Iy += func(x+hx,y+hy);
				y += hy;
			}

			I += Iy * hy / 6;
		}

		x += hx;
	}

	return I * hx / 6;	
}

// =========================================================================

Simpson3D::Simpson3D ( const vector_tuple_did &limits )
: x_min(std::get<0>(limits[0])), x_terms(std::get<1>(limits[0])), x_max(std::get<2>(limits[0])),
  y_min(std::get<0>(limits[1])), y_terms(std::get<1>(limits[1])), y_max(std::get<2>(limits[1])),
  z_min(std::get<0>(limits[2])), z_terms(std::get<1>(limits[2])), z_max(std::get<2>(limits[2])) 
{ 
	if (limits.size() != 3) 
		throw std::invalid_argument("Simpson3D: Only 3dim is implemented!");
	if ( (x_min >= x_max) || (y_min >= y_max) || (z_min >= z_max) )
		throw std::invalid_argument("Low bound of integral is bigger then upper.");
}

double Simpson3D::value (const std::function<double(double,double,double)> &func) const
{
	double hx = (this->x_max - this->x_min) / this->x_terms;
	double hy = (this->y_max - this->y_min) / this->y_terms;
	double hz = (this->z_max - this->z_min) / this->z_terms;

	double total_volume = (this->x_max-this->x_min) 
						* (this->y_max-this->y_min) 
						* (this->z_max-this->z_min);

	double total_terms = this->x_terms * this->y_terms * this->z_terms;
	double I = 0;

	std::vector<std::vector<double>> 
	arg_grid(3, std::vector<double>
			(3, 0.0));

	arg_grid[0][0] = this->x_min;							// x_a
	arg_grid[0][2] = this->x_min + hx;						// x_ab
	arg_grid[0][1] = (arg_grid[0][0] + arg_grid[0][2]) / 2;	// x_b

	arg_grid[1][0] = this->y_min;							// y_a
	arg_grid[1][2] = this->y_min + hy;						// y_ab
	arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2; // y b

	arg_grid[2][0] = this->z_min;							// z_a
	arg_grid[2][2] = this->z_min + hz;						// z_ab
	arg_grid[2][1] = (arg_grid[2][0] + arg_grid[2][2]) / 2; // z_b

	// init setup for fun_grid
	std::vector<std::vector<std::vector<double>>>
	fun_grid(3, std::vector<std::vector<double>>
			(3, std::vector<double>
			(3, 0.0)));

	for (std::size_t x = 0; x < 3; x++) {
		for (std::size_t y = 0; y < 3; y++) {
			for (std::size_t z = 0; z < 3; z++) {
				fun_grid[x][y][z] = func(arg_grid[0][x],
										 arg_grid[1][y],
										 arg_grid[2][z]);
			}
		}
	}

	// init setup for c_grid
	std::vector<std::vector<std::vector<double>>>
	c_grid  (3, std::vector<std::vector<double>>
			(3, std::vector<double>
			(3, 1.0)));

	for (std::size_t x = 0; x < 3; x++) {
		for (std::size_t y = 0; y < 3; y++) {
			for (std::size_t z = 0; z < 3; z++) {
				if (x == 1) c_grid[x][y][z] *= 4.0;
				if (y == 1) c_grid[x][y][z] *= 4.0;
				if (z == 1) c_grid[x][y][z] *= 4.0;
			}
		}
	}

	while (arg_grid[0][2] <= this->x_max) {

		arg_grid[1][0] = this->y_min;
		arg_grid[1][2] = this->y_min + hy;
		arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2;
		
		while (arg_grid[1][2] <= this->y_max) {

			arg_grid[2][0] = this->z_min;
			arg_grid[2][2] = this->z_min + hz;
			arg_grid[2][1] = (arg_grid[2][0] + arg_grid[2][2]) / 2;
			
			while (arg_grid[2][2] <= this->z_max) {

				arg_grid[2][0] += hz; // UPD: z_from
				arg_grid[2][1] += hz; // UPD: z_mid
				arg_grid[2][2] += hz; // UPD: z_to

				for (std::size_t X = 0; X < 3; X++) {
					for (std::size_t Y = 0; Y < 3; Y++) {
						for (std::size_t Z = 0; Z < 3; Z++) {
							// summation of 27 terms
							I += c_grid[X][Y][Z] * fun_grid[X][Y][Z];
							// update fun_grid matrix
							// TODO: not optimum
							fun_grid[X][Y][Z] = func(arg_grid[0][X],
													 arg_grid[1][Y],
													 arg_grid[2][Z]);
						}
					}
				}

			}
			arg_grid[1][0] += hy; // UPD: y_from
			arg_grid[1][1] += hy; // UPD: y_mid
			arg_grid[1][2] += hy; // UPD: y_to
		}
		arg_grid[0][0] += hx; // UPD: x_from
		arg_grid[0][1] += hx; // UPD: x_mid
		arg_grid[0][2] += hx; // UPD: x_to
	}

	return total_volume / total_terms / 216.0 * I;
}

// =========================================================================

Simpson2D::Simpson2D ( const vector_tuple_did &limits )
: x_min(std::get<0>(limits[0])), x_terms(std::get<1>(limits[0])), x_max(std::get<2>(limits[0])),
  y_min(std::get<0>(limits[1])), y_terms(std::get<1>(limits[1])), y_max(std::get<2>(limits[1])) 
{ 
	if (limits.size() != 2) 
		throw std::invalid_argument("Simpson2D: Only 3dim is implemented!");
	if ( (x_min >= x_max) || (y_min >= y_max) )
		throw std::invalid_argument("Low bound of integral is bigger then upper.");
}

double Simpson2D::value (const std::function<double(double,double)> &func) const
{
	double hx = (this->x_max - this->x_min) / this->x_terms;
	double hy = (this->y_max - this->y_min) / this->y_terms;

	double total_volume = (this->x_max-this->x_min) 
						* (this->y_max-this->y_min);

	double total_terms = this->x_terms * this->y_terms;
	double I = 0;

	std::vector<std::vector<double>> 
	arg_grid(2, std::vector<double>
			(3, 0.0));

	arg_grid[0][0] = this->x_min;							// x_a
	arg_grid[0][2] = this->x_min + hx;						// x_ab
	arg_grid[0][1] = (arg_grid[0][0] + arg_grid[0][2]) / 2;	// x_b

	arg_grid[1][0] = this->y_min;							// y_a
	arg_grid[1][2] = this->y_min + hy;						// y_ab
	arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2; // y_b

	// init setup for fun_grid
	std::vector<std::vector<double>>
	fun_grid(3, std::vector<double>
			(3, 0.0));

	for (std::size_t x = 0; x < 3; x++) {
		for (std::size_t y = 0; y < 3; y++) {
			fun_grid[x][y] = func(arg_grid[0][x], arg_grid[1][y]);
		}
	}

	// init setup for c_grid
	std::vector<std::vector<double>>
	c_grid  (3, std::vector<double>
			(3, 1.0));

	for (std::size_t x = 0; x < 3; x++) {
		for (std::size_t y = 0; y < 3; y++) {
			if (x == 1) c_grid[x][y] *= 4.0;
			if (y == 1) c_grid[x][y] *= 4.0;
		}
	}

	while (arg_grid[0][2] <= this->x_max) {

		arg_grid[1][0] = this->y_min;
		arg_grid[1][2] = this->y_min + hy;
		arg_grid[1][1] = (arg_grid[1][0] + arg_grid[1][2]) / 2;
		
		while (arg_grid[1][2] <= this->y_max) {

			arg_grid[1][0] += hy; // UPD: y_from
			arg_grid[1][1] += hy; // UPD: y_mid
			arg_grid[1][2] += hy; // UPD: y_to

			for (std::size_t X = 0; X < 3; X++) {
				for (std::size_t Y = 0; Y < 3; Y++) {
					I += c_grid[X][Y] * fun_grid[X][Y];
					fun_grid[X][Y] = func(arg_grid[0][X], arg_grid[1][Y]);
				}
			}

		}

		arg_grid[0][0] += hx; // UPD: x_from
		arg_grid[0][1] += hx; // UPD: x_mid
		arg_grid[0][2] += hx; // UPD: x_to
	}

	return total_volume / total_terms / 36.0 * I;
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