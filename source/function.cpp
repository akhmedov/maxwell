//
//  function.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 18.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "function.hpp"

double Function::rect (double x, double duration)
{
	return (x >= 0 && x <= std::abs(duration)) ? 1 : 0;
}

double Function::sin (double x, double duration)
{
	duration = std::abs(duration);
	if (x >= 0 && x <= duration)
		return std::sin(M_PI * x / duration);
	return 0;
}

double Function::sinc (double x, double duration, std::size_t cycles)
{
	duration = std::abs(duration);
	if (x >= 0 && x <= duration) {
			x -= duration / 2;
			x *= 2 * cycles * M_PI / duration;
			return (x) ? std::sin(x)/x : 1;
	} return 0;
}

double Function::gauss (double x, double duration)
{
	double mu = std::abs(duration) / 2;
	double sigma = 2 * std::pow(std::abs(duration)/7,2);;
	double pow = std::pow(x-mu,2) / sigma;
	return std::exp(-pow);
}

double Function::gauss_perp (double x, double duration, std::size_t order)
{
	double a = duration / 2.0;
	double b = std::pow(duration / 7,2);
	double c1 = (x - a) / b;
	double c2 = 1.0 / b;

	switch (order) {
		case 0: {
			return Function::gauss(x, duration);
		} case 1: {
			return - c1 * Function::gauss(x, duration);
		} default: {
			double parent = Function::gauss_perp(x, duration, order-1);
			double grand = Function::gauss_perp(x, duration, order-2);
			double drivative = - (order-1.0) * c2 * grand - c1 * parent;
			return drivative;
		}
	}
}

double Function::gauss_perp_normed (double x, double duration, std::size_t order)
{
	duration -= order*duration/20;
	x -= order*duration/40;
	double norm1 = tgamma(order/2+1) / tgamma(order+1);
	double norm2 = std::pow(duration/7,order);
	return norm1 * norm2 * Function::gauss_perp(x, duration, order);
}

double Function::sigmoid (double x, double duration)
{	
	if (x >= 0 && x <= duration) {
		x -= duration / 2;
		x *= 12 / duration;
		return 1 / (1 + std::exp(-x));
	}
	
	if (x > duration) return 1;
	return 0;
}

double Function::smoozed_rect (double x, double duration, double sensivity)
{
	duration -= sensivity;
	return Function::sigmoid (x, sensivity) - Function::sigmoid (x - duration, sensivity);	
}
