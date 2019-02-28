//
//  vector_to_matrix.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 27.02.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
#include "gnu_plot.hpp"

#include <iomanip>
#include <iostream>
using namespace std;

struct Test {
    static bool vector_to_matrix ();
};

bool Test::vector_to_matrix ()
{
	std::vector<std::vector<double>> ARRAY, MATRIX, matrix;

	auto equals = [] (const std::vector<std::vector<double>> &a, 
					  const std::vector<std::vector<double>> &b) {
		for (std::size_t i = 0; i < a.size(); i++) {
			auto inv_i = a.size() - i - 1;
			for (std::size_t j = 0; j < a[0].size(); j++)
				if (a[inv_i][j] != b[i][j]) return false;
		} return true;
	};

	ARRAY = {
		{-1, 1, 3}, { 0, 1, 0}, { 1, 1, 3},
		{-1, 0, 3}, { 0, 0, 3}, { 1, 0, 0},
		{-1,-1, 0}, { 0,-1, 0}, { 1,-1, 3}
	};

	MATRIX = {
		{3, 0, 3},
		{3, 3, 0},
		{0, 0, 3}
	};

	matrix = GnuPlot::datagrid_from (ARRAY, 0, 1);
	if (!equals(MATRIX,matrix)) return false;

	ARRAY = {
		{-1, 1, 3}, { 0, 1, 0}, { 1, 1, 3},
		{-1, 0, 3}, { 0, 0, 3}, { 1, 0, 0}
	};

	MATRIX = {
		{3, 0, 3},
		{3, 3, 0}
	};

	matrix = GnuPlot::datagrid_from (ARRAY, 0, 1);
	if (!equals(MATRIX,matrix)) return false;

	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::vector_to_matrix()" << left;
    cout << (Test::vector_to_matrix() ? "PASSED" : "FAILED") << endl;
}
