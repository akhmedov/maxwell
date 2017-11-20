//
//  nonlinear_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef nonlinear_field_hpp
#define nonlinear_field_hpp

#include "abstract_field.hpp"
#include "linear_field.hpp"
#include "nonlinear_medium.hpp"

#include <complex>

struct NonlinearField : public AbstractField {

	NonlinearField (LinearField* field, NonlinearMedium* medium);

	virtual double electric_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double electric_z (double ct, double rho, double phi, double z) const = 0;

	virtual double magnetic_rho (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_phi (double ct, double rho, double phi, double z) const = 0;
	virtual double magnetic_z (double ct, double rho, double phi, double z) const = 0;

protected:
	LinearField* field;
	NonlinearMedium* nl_medium;
};

#endif /* nonlinaer_field_hpp */
