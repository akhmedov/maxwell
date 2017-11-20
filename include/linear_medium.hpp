//
//  linear_medium.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 26.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef linear_medium_hpp
#define linear_medium_hpp

#include <cmath>

struct LinearMedium {
	virtual double relative_permittivity (double ct, double z) const = 0;
	virtual double relative_permeability (double ct, double z) const = 0;
	double signal_max_velocity (double ct, double z) const;
	
	// TODO: use GMP lib double for next constants
	constexpr static const double C = 299792458;
	constexpr static const double C2 = C * C;
	constexpr static const double EPS0 = 10e7 / (4 * M_PI * C2);
	constexpr static const double MU0 = 4 * M_PI * 10e-7;
};

#endif /* linear_medium_hpp */
