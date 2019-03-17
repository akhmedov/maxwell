//
//  snr_dataset.cpp
//  example.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "dataset.hpp"

static const int 	RADIX  = 3;
static const float 	DCYCLE = 0.5;
static const float 	NPOWER = 10;

#include <vector>
#include <iostream>
using namespace std;

int main ()
{
	Dataset* ds = new Dataset(RADIX, DCYCLE, NPOWER);
    return 0;
}
