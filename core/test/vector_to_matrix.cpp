//
//  vector_to_matrix.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 27.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "script_manager.hpp"
using namespace std;

struct Test : public ScriptManager {
    static bool vector_to_matrix ();
};

bool Test::vector_to_matrix ()
{
	std::vector<std::pair<double,double>> XY;
	std::vector<double> Z;
	std::vector<std::vector<double>> MATRIX, matrix;
	double xmin, xstep, xmax, ymin, ystep, ymax;

	auto equals = [] (const std::vector<std::vector<double>> &a, 
					  const std::vector<std::vector<double>> &b) {
		for (std::size_t i = 0; i < a.size(); i++) {
			auto inv_i = a.size() - i - 1;
			for (std::size_t j = 0; j < a[0].size(); j++)
				if (a[inv_i][j] != b[i][j]) return false;
		} return true;
	};

	XY = {{-1,1}, {0,1}, {1,1}, {-1,0}, {0,0}, {1,0}, {-1,-1}, {0,-1}, {1,-1}};
	Z = {3, 0, 3, 3, 3, 0, 0, 0, 3};

	MATRIX = {
		{3, 0, 3},
		{3, 3, 0},
		{0, 0, 3}
	};

	matrix = ScriptManager::datagrid_from(XY, Z, xmin, xstep, xmax, ymin, ystep, ymax);
	if (!equals(MATRIX,matrix)) return false;

	XY = {{-1,1}, {0,1}, {1,1}, {-1,0}, {0,0}, {1,0}};
	Z = {3, 0, 3, 3, 3, 0};

	MATRIX = {
		{3, 0, 3},
		{3, 3, 0}
	};

	matrix = ScriptManager::datagrid_from(XY, Z, xmin, xstep, xmax, ymin, ystep, ymax);
	if (!equals(MATRIX,matrix)) return false;

	return true;
}

int main ()
{
	return !Test::vector_to_matrix();
}
