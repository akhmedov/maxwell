//
//  nonlinear_medium.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "nonlinear_medium.hpp"

double NonlinearMedium::relative_permittivity (double ct, double z) const
{
	return this->relative_permittivity(ct, z, 1);
}

double NonlinearMedium::relative_permeability (double ct, double z) const
{
	return this->relative_permeability(ct, z, 1);
}
