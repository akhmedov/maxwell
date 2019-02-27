//
//  nonlinear_medium.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef nonlinear_medium_hpp
#define nonlinear_medium_hpp

#include "linear_medium.hpp"

#include <vector>

struct NonlinearMedium : public LinearMedium {
	virtual double conductivity (double ct, double z) const = 0;
	double relative_permittivity (double ct, double z) const;
	double relative_permeability (double ct, double z) const;

	/* TODO: implement next form of call: 		   *
	 * anizotropic.relative_permittivity<3>(0, 20) */
	virtual double relative_permittivity (double ct, double z, std::size_t term) const = 0;
	virtual double relative_permeability (double ct, double z, std::size_t term) const = 0;

};

#endif /* nonlinear_medium_hpp */
