//
//  I1_time_partder.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 27.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "phys_math.hpp"
#include "integral.hpp"

#include "updisk_meandr.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

using namespace std;

struct Test : KerrAmendment {
    static bool I1_time_partder ();
};

bool Test::I1_time_partder ()
{
	for (double rho = 0; rho <= 1.3; rho += 0.25) {
		for (double R = 0.1; R <= 1.3; R += 0.2) {
			for (double ct = 0.1; ct <= 1.4; ct += 0.3) {
				for (double z = 1e-2; z <= ct - 0.05; z += 0.3) {

					auto I1 = [rho, R, z] (double ct_perp) {
						double vt_z = std::sqrt(ct_perp * ct_perp - z * z);
						return TransientResponse::int_bessel_011(vt_z,rho,R);
					};
					
					double vt_z = ct * ct - z * z;
					double anal = ct * KerrAmendment::int_bessel_011_perp(vt_z,rho,R);
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
	return !Test::I1_time_partder();
}
