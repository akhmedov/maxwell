//
//  anal_numerical_I1.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 27.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include <iomanip>
#include <iostream>
using namespace std;

bool anal_numerical_I1 ()
{
    return false;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::anal_numerical_I1()" << left;
    cout << (anal_numerical_I1() ? "PASSED" : "FAILED") << endl;
}
