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

#include <random>

struct Noise {
	virtual double value (double filed, double ct, double rho, double phi, double z) = 0;
};

struct AdditiveWhite : public Noise {
	AdditiveWhite (double magnutude, double diviation);
	double value (double filed, double ct, double rho, double phi, double z);
private:
	std::uniform_real_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double diviation;
};

struct AdditiveWhiteGaussian : public Noise {
	AdditiveWhiteGaussian (double magnutude, double sigma);
	double value (double filed, double ct, double rho, double phi, double z);
private:
	std::normal_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double sigma;
};

struct MultiplicWhiteGaussian : public Noise {
	MultiplicWhiteGaussian (double magnutude, double sigma);
	double value (double filed, double ct, double rho, double phi, double z);
private:
	std::normal_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double sigma;
};

#endif /* noise_hpp */
