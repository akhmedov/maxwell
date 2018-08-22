//
//  noise.hpp
//  Evolution
//
//  Created by Rolan Akhmedov on 31.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#ifndef noise_hpp
#define noise_hpp

#include <abstract_field.hpp>

#include <random>

struct Noise {
	virtual double value (double ct, double rho, double phi, double z, double filed = 0) = 0;
	virtual double power (double rho, double phi, double z, std::size_t samples = 1e3);
};

struct AdditiveWhite : public Noise {
	AdditiveWhite (double magnutude, double diviation);
	double value (double ct, double rho, double phi, double z, double filed = 0);
private:
	std::uniform_real_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double diviation;
};

struct AdditiveWhiteGaussian : public Noise {
	AdditiveWhiteGaussian (double magnutude, double sigma);
	double value (double ct, double rho, double phi, double z, double filed = 0);
private:
	std::normal_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double sigma;
};

struct MultiplicWhiteGaussian : public Noise {
	MultiplicWhiteGaussian (double magnutude, double sigma);
	double value (double ct, double rho, double phi, double z, double filed = 1);
private:
	std::normal_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double sigma;
};

struct NoiseField : public AbstractField {

	NoiseField (double magnutude, double sigma);
	double electric_rho (double ct, double rho, double phi, double z) const;
	double electric_phi (double ct, double rho, double phi, double z) const;
	double electric_z (double ct, double rho, double phi, double z) const;
	double magnetic_rho (double ct, double rho, double phi, double z) const;
	double magnetic_phi (double ct, double rho, double phi, double z) const;
	double magnetic_z (double ct, double rho, double phi, double z) const;

private:

	Noise* radial;
	Noise* athimus;
	Noise* distance;
};

#endif /* noise_hpp */
