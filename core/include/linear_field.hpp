//
//  linear_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef linear_field_hpp
#define linear_field_hpp

#include "abstract_field.hpp"
#include "linear_source.hpp"
#include "linear_medium.hpp"

struct LinearField : AbstractField {

	LinearField (LinearCurrent* linear_source, LinearMedium* linear_medium);
	
	virtual double electric_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_z (double ct, double rho, double phi, double z) const = 0;

	virtual double magnetic_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_z (double ct, double rho, double phi, double z) const = 0;
	
protected:
	LinearCurrent* source;
	LinearMedium* medium;
};

#endif /* linaer_field_hpp */
