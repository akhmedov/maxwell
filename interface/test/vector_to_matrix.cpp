//
//  vector_to_matrix.cpp
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

bool UnitTest::hardcode_dataset ()
{
	double tau = 0.5;
	double duty_cycle = 0.5;
	double vt_step = 0.01;

	std::pair<double,double> rho = std::make_pair(0,7);
	std::pair<double,double> phi = std::make_pair(0,90);
	std::pair<double,double> z = std::make_pair(0,20);

	std::vector<std::function<double(double)>> domain = {
		[tau] (double vt) { return   Function::sinc  (vt,tau); },
		[tau] (double vt) { return   Function::gauss (vt,tau); }
	};

	serial::dataset ds = serial::dataset(domain, tau, NULL, duty_cycle, vt_step);
	
	ds.set_char(0, 0, 2, 1);
	ds.set_char(0, 0, 2, 2);
	ds.set_char(0, 0, 2, 2);
	ds.set_char(0, 0, 2, 1);

	json js = serial::json_from(ds,rho,phi,z);
	serial::serialize("dataset.json",js);
	return true;
}


bool UnitTest::random_dataset ()
{
	double radix = 3;
	double sigma = 5.0;
	double pulses = 60;
	std::string file_name = "dataset.json";
	std::pair<double,double> rho = std::make_pair(0,4);
	std::pair<double,double> phi = std::make_pair(0,20);
	std::pair<double,double> z = std::make_pair(4,6);
	serial::randomized_sequental (pulses, radix, sigma, file_name, rho, phi, z);
	return true;
}

bool UnitTest::snr_dataset ()
{
	double radix = 3;
	double snr = 40;
	double pulses = 1e3;
	std::string file_name = "train.json";
	double rho = 0;
	double phi = 0;
	double z = 10;
	serial::same_snr(pulses, radix, snr, file_name, rho, phi, z);
	return true;
}

int main ()
{
	cout << left << setfill('.') << setw(70);
	cout << "Test::Math::vector_to_matrix()" << left;
    cout << (Test::vector_to_matrix() ? "PASSED" : "FAILED") << endl;
}
