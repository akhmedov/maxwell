//
//  noise.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 31.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "noise.hpp"


double Noise::power (double rho, double phi, double z, std::size_t samples)
{
	double sum = 0;
	for (double vt = 0; vt <= samples * 0.01; vt += 0.01) {
		double noise = this->value(vt,rho,phi,z);
		sum += noise * noise;
	}
	return sum / samples;
}

//=============================================================================


AdditiveWhite::AdditiveWhite (double avarage_magnutude, double diviation_magnutude)
{ 
	this->magnitude = avarage_magnutude;
	this->diviation = diviation_magnutude;

	double from = avarage_magnutude - diviation_magnutude/2;
	double to = avarage_magnutude + diviation_magnutude/2;
	this->distribution = std::uniform_real_distribution<double>(from,to);
	this->generator = std::mt19937_64(this->device());
}

double AdditiveWhite::value (double ct, double rho, double phi, double z, double filed)
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

double AdditiveWhiteGaussian::value (double ct, double rho, double phi, double z, double filed)
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

double MultiplicWhiteGaussian::value (double ct, double rho, double phi, double z, double filed)
{
	UNUSED(ct); UNUSED(rho); UNUSED(phi); UNUSED(z);
	return filed * this->distribution(this->generator);
}

//============================================================================

NoiseField::NoiseField (double magnutude, double sigma)
{
	#ifdef UNIFORM_NOICE
		this->radial = new AdditiveWhite(magnutude, sigma);
		this->athimus = new AdditiveWhite(magnutude, sigma);
		this->distance = new AdditiveWhite(magnutude, sigma);		
	#else
		this->radial = new AdditiveWhiteGaussian(magnutude, sigma);
		this->athimus = new AdditiveWhiteGaussian(magnutude, sigma);
		this->distance = new AdditiveWhiteGaussian(magnutude, sigma);
	#endif /* UNIFORM_NOICE */
}

double NoiseField::electric_rho (double ct, double rho, double phi, double z) const
{
	return this->radial->value(ct, rho, phi, z);
}

double NoiseField::electric_phi (double ct, double rho, double phi, double z) const
{
	return this->athimus->value(ct, rho, phi, z);
}

double NoiseField::electric_z (double ct, double rho, double phi, double z) const
{
	return this->distance->value(ct, rho, phi, z);
}

double NoiseField::magnetic_rho (double ct, double rho, double phi, double z) const
{
	return this->radial->value(ct, rho, phi, z);
}

double NoiseField::magnetic_phi (double ct, double rho, double phi, double z) const
{
	return this->athimus->value(ct, rho, phi, z);
}

double NoiseField::magnetic_z (double ct, double rho, double phi, double z) const
{
	return this->distance->value(ct, rho, phi, z);
}
