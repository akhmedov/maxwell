//
//  linear_source.hpp
//  Evolution
//
//  Created by Rolan Akhmedov on 26.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef linear_source_hpp
#define linear_source_hpp

#include <vector>
#include <cmath>
#include <vector>

struct LinearSource { };

struct LinearCurrent : public LinearSource {
	virtual double rho (double ct, double rho, double phi, double z) const = 0;
	virtual double phi (double ct, double rho, double phi, double z) const = 0;
	virtual double z (double ct, double rho, double phi, double z) const = 0;

	double x (double ct, double rho, double phi, double z) const;
	double y (double ct, double rho, double phi, double z) const;
	std::vector<double> cylindric (double ct, double rho, double phi, double z) const;
	std::vector<double> cartesian (double ct, double rho, double phi, double z) const;
};

struct LinearCharge : public LinearSource {
	virtual double value (double ct, double rho, double phi, double z) const = 0;
};

#endif /* linear_source_hpp */
