//
//  electrodynamics.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"
#include "abstract_field.hpp"

struct Electrodynamics : public AbstractField {
	double directivity () const;
	double energy_e (double rho, double phi, double z) const;
	double energy_eh (double rho, double phi, double z) const;
};
