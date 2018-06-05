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
	if (x >= 0 && x <= duration) {
		double mu = std::abs(duration) / 2;
		double sigma = 2 * std::pow(std::abs(duration)/7,2);;
		double pow = std::pow(x-mu,2) / sigma;
		return std::exp(-pow); // / std::sqrt(M_PI * sigma);
	} return 0;
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
