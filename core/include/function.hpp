//
//  function.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 18.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef function_hpp
#define function_hpp

#include "phys_math.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace Function {
	double sin (double x, double duration);
	double rect (double x, double duration);
	double sinc (double x, double duration, std::size_t cycles = 10);
	double gauss (double x, double duration);
	double gauss_perp (double x, double duration, std::size_t order);
	double gauss_perp_normed (double x, double duration, std::size_t order);
	double sigmoid (double x, double duration);
	double smoozed_rect (double x, double duration, double sensivity);
	// double gauss_perp (double x, double duration, std::size_t order);
};

#endif /* function_hpp */
