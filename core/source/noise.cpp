//
//  noise.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 31.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "noise.hpp"


double Noise::power (std::size_t samples)
{
	double sum = 0;
	for (double vt = 0; vt <= samples * 0.01; vt += 0.01) {
		double noise = this->value();
		sum += noise * noise;
	}
	return sum / samples;
}

//=============================================================================


WhiteUniform::WhiteUniform (double avarage_magnutude, double diviation_magnutude)
{ 
	this->magnitude = avarage_magnutude;
	this->diviation = diviation_magnutude;

	double from = avarage_magnutude - diviation_magnutude/2;
	double to = avarage_magnutude + diviation_magnutude/2;
	this->distribution = std::uniform_real_distribution<double>(from,to);
	this->generator = std::mt19937_64(this->device());
}

double WhiteUniform::value ()
{
	return this->distribution(this->generator);
}

//=============================================================================

WhiteGaussian::WhiteGaussian (double avarage_magnutude, double sigma_magnutude)
{
	this->magnitude = avarage_magnutude;
	this->sigma = sigma_magnutude;

	this->distribution = std::normal_distribution<double>(avarage_magnutude, sigma_magnutude);
	this->generator = std::mt19937_64(this->device());
}

double WhiteGaussian::value ()
{
	return this->distribution(this->generator);
}
