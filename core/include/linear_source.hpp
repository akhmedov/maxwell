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

struct LinearSource { 
	LinearSource (double tau = NAN);
	double get_duration () const;
	virtual double time_shape (double vt) const = 0;
protected:
	const double duration;
};

struct LinearCurrent : public LinearSource {

	LinearCurrent (double tau = NAN);

	virtual double rho (double rho, double phi, double z) const = 0;
	virtual double phi (double rho, double phi, double z) const = 0;
	virtual double z   (double rho, double phi, double z) const = 0;

	double x (double rho, double phi, double z) const;
	double y (double rho, double phi, double z) const;
	std::vector<double> cylindric (double rho, double phi, double z) const;
	std::vector<double> cartesian (double rho, double phi, double z) const;
};

struct LinearCharge : public LinearSource {
	LinearCharge (double tau = NAN);
	virtual double value (double rho, double phi, double z) const = 0;
};

#endif /* linear_source_hpp */
