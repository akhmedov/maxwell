//
//  noise.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 31.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "noise.hpp"

/* std::uniform_real_distribution<double> distribution;
std::mt19937_64 generator;
double magnitude;
double value; */

AdditiveWhite::AdditiveWhite (double avarage_magnutude, double diviation_magnutude)
{ 
	this->magnitude = avarage_magnutude;
	this->diviation = diviation_magnutude;

	double from = avarage_magnutude - diviation_magnutude/2;
	double to = avarage_magnutude + diviation_magnutude/2;
	this->distribution = std::uniform_real_distribution<double>(from,to);
	this->generator = std::mt19937_64(this->device());
}

double AdditiveWhite::value (double filed, double ct, double rho, double phi, double z)
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	return filed + this->distribution(this->generator);
}

//=============================================================================

AdditiveWhiteGaussian::AdditiveWhiteGaussian (double avarage_magnutude, double sigma_magnutude)
{
	this->magnitude = avarage_magnutude;
	this->sigma = sigma_magnutude;

	this->distribution = std::normal_distribution<double>(avarage_magnutude, sigma_magnutude);
	this->generator = std::mt19937_64(this->device());
}

double AdditiveWhiteGaussian::value (double filed, double ct, double rho, double phi, double z)
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	return filed + this->distribution(this->generator);
}

//=============================================================================

MultiplicWhiteGaussian::MultiplicWhiteGaussian (double avarage_magnutude, double sigma_magnutude)
{
	this->magnitude = avarage_magnutude;
	this->sigma = sigma_magnutude;

	this->distribution = std::normal_distribution<double>(avarage_magnutude, sigma_magnutude);
	this->generator = std::mt19937_64(this->device());
}

double MultiplicWhiteGaussian::value (double filed, double ct, double rho, double phi, double z)
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	return filed * this->distribution(this->generator);
}
