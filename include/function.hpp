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
#include <stdexcept>

namespace Function {
	double sin (double x, double duration);
	double rect (double x, double duration);
	double sinc (double x, double duration, std::size_t cycles = 10);
	double gauss (double x, double duration);
	double sigmoid (double x, double duration);
};

#endif /* function_hpp */
