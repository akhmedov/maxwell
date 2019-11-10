//
//  fft.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 10.11.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"
#include "script_manager.hpp"

#include <vector>
#include <complex>
#include <algorithm>

using namespace std;

struct Test : Math {
    static bool harmonic_function ();
};

bool Test::harmonic_function ()
{
    double w = 3;
    auto fnc = [w] (double t) { return sin(w*t); };
    
    vector<complex<double>> fncVal;
    for (double arg = 0; arg < 2*M_PI/w; arg += 2*M_PI/10/w)
        fncVal.push_back(fnc(arg));

    vector<complex<double>> imgVal = Math::fft(fncVal);
    for_each(imgVal.begin(), imgVal.end(), [](complex<double> &n){ n /= 10; });
    vector<complex<double>> newFncVal = Math::inv_fft(imgVal);

    for (size_t i = 0; i < newFncVal.size(); i++) {
        if (abs(fncVal[i]) < 1e-10) continue;
        double error = 100 * abs(newFncVal[i] - fncVal[i]) / abs(fncVal[i]);
        if (error > 1) return false;
    }

    return true;
}

int main ()
{
	return !Test::harmonic_function();
}
