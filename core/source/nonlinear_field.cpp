//
//  nonlinear_field.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 22.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "nonlinear_field.hpp"

NonlinearField::NonlinearField (LinearField* linear_field, NonlinearMedium* medium)
: field(linear_field), nl_medium(medium) { }
