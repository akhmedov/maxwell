//
//  linear_medium.cpp
//  Evolution
//
//  Created by Rolan Akhmedov on 26.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "linear_medium.hpp"

double LinearMedium::signal_max_velocity (double ct, double z) const
{
	const double mu_r = this->relative_permeability(ct, z);
	const double eps_r = this->relative_permittivity(ct, z);
	const double sqrt_mu_eps = std::sqrt(mu_r * eps_r);
	return LinearMedium::C / sqrt_mu_eps;
}

