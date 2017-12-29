//
//  linear_field.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "linear_field.hpp"

LinearField::LinearField (LinearCurrent* linear_source, LinearMedium* linear_medium)
: source(linear_source), medium(linear_medium) { }

