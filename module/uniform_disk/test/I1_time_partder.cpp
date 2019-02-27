//
//  I1_time_partder.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 27.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "integral.hpp"

#include <iomanip>
#include <iostream>
using namespace std;

struct Test {
    static bool I1_time_partder ();
};

bool Test::I1_time_partder ()
{
	for (double rho = 0; rho <= 1.3; rho += 0.25) {
		for (double R = 0.1; R <= 1.3; R += 0.2) {
			for (double ct = 0.1; ct <= 1.4; ct += 0.3) {
				for (double z = 1e-2; z <= ct - 0.05; z += 0.3) {

					auto I1 = [rho, R, z] (double ct) {
						double vt_z = std::sqrt(ct * ct - z * z);
						return MissileField::int_bessel_011(vt_z,rho,R);
					};
					
					double anal = KerrAmendment::int_bessel_011_perp(ct,z,rho,R);
					auto num3 = Math::derivat3(I1, ct);
					auto num4 = Math::derivat4(I1, ct);

					double error3 = anal ? std::abs(anal - num3) / anal : 0;
					double error4 = anal ? std::abs(anal - num4) / anal : 0;
					if (100*error3 > 10) return false;
					if (100*error4 > 5) return false;
				}
			}
		}
	}

	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::I1_time_partder()" << left;
    cout << (Test::I1_time_partder() ? "PASSED" : "FAILED") << endl;
}
