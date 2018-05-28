//
//  electrodynamics.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef electrodynamics_hpp
#define electrodynamics_hpp

#include "integral.hpp"
#include "abstract_field.hpp"

struct Electrodynamics : public AbstractField {
	double energy_cart (double x, double y, double z) const;
	double energy (double rho, double phi, double z) const;
};

#endif /* electrodynamics_hpp */
