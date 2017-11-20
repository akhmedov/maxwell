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

double Simpson::value (double from, double to, std::function<double(double)> func)
{
    if (from >= to) throw std::invalid_argument("Low bound of integral is bigger then upper.");

	double h = (to - from) / this->quadr_terms;
	double I = 0;
	for (std::size_t i = 0; i < this->quadr_terms; i++) {
		double a = from + h * i;
		double b = from + h * (i+1);
		I += func(a) + 4 * func((a + b)/2) + func(b);
	}

	I = I * (to - from) / 6;
	return I / this->quadr_terms;
}

// =========================================================================

GaussLaguerre::GaussLaguerre ()
: quadr_terms(polynom_data.size()) { }

GaussLaguerre::GaussLaguerre (std::size_t terms)
: quadr_terms(polynom_data.size()) {
	if (terms > this->polynom_data.size()) 
		throw std::logic_error("Max terms number extended!");
}

double GaussLaguerre::value (std::function<double(double)> func)
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
