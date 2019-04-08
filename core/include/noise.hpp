//
//  noise.hpp
//  Evolution
//
//  Created by Rolan Akhmedov on 31.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include <random>

struct Noise {
	virtual double value () = 0;
	virtual double power (std::size_t samples = 1e3);
};

struct WhiteUniform : public Noise {
	WhiteUniform (double magnutude, double diviation);
	double value ();
private:
	std::uniform_real_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double diviation;
};

struct WhiteGaussian : public Noise {
	WhiteGaussian (double magnutude, double sigma);
	double value ();
private:
	std::normal_distribution<double> distribution;
	std::random_device device;
	std::mt19937_64 generator;
	double magnitude;
	double sigma;
};
