//
//  main.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include <cfloat> // size_t
#include <iostream> // cout cerr
#include <string> // compare substr
#include <vector> // push_back

#include "ready_model.hpp"
#include "phys_math.hpp"
#include "config.hpp"

using namespace std;

vector<float> cstring_to_string (size_t from, size_t to, char* argv[]);
void print_configuration ();

int main(int argc, char* argv[])
{

	// find entity of contig option
	string posix_path_conf = "maxwell.conf";
	for (int j = 1; j < argc; j++) {
		string iter = string(argv[j]);
		if (iter == "--conf") {
			posix_path_conf = string(argv[++j]);
			break;
		}
	}

	int i = 1;
	string item = string(argv[i]);

	if ( item == "--help" ) {
		cout << "--help \t\t\t show this menu" << endl;
		cout << "--version \t\t build information" << endl;
		cout << "--conf <path> \t\t maxwell configuration file" << endl;
		cout << "--model <num> \t\t plot model number" << endl;
		cout << "\t num 1: \t linear Ex(ct)" << endl;
		cout << "\t\t\t --arg <rho> <phi>" << endl;
		cout << "\t num 2: \t linear Hy(ct)" << endl;
		cout << "\t\t\t --arg <z>" << endl;
		cout << "\t num 3: \t linear Hy(ct, z)" << endl;
		cout << "\t\t\t --arg <terms_num>" << endl;
		cout << "\t num 4: \t linear Ex(ct, z)" << endl;
		cout << "\t num 5: \t linear Ex(ct, phi)" << endl;
		cout << "\t\t\t --arg <rho>" << endl;
		cout << "\t num 6: \t linear Ex(ct, rho)" << endl;
		cout << "\t num 7: \t missile effect length on oZ" << endl;
		cout << "\t num 8: \t magnetic static magnitude from oZ" << endl;
		cout << "\t num 9: \t magnetic static magnitude with different terms number" << endl;
		cout << "\t num 10: \t numerical integration compare for I1" << endl;
		cout << "\t num 11: \t numerical integration compare for I2" << endl;
		cout << "--arg <val> ... <val> \t arguments for plot model" << endl;
		return 0;
	}

	else if ( item == "--version" ) {
		cout << "Maxwell 0.9 alpha" << endl;
		cout << "Copyright (C) 2017 Rolan Akhmedov (Ukraine)" << endl;
		cout << "Contact e-mail: rolan.kharkiv@gmail.com" << endl;
		return 0;
	}

	else if ( item == "--model" ) {
		cout << "Configuration file: " << posix_path_conf << endl;
		Config::read(posix_path_conf);
		if ( Config::display_params() ) print_configuration();
		int model_number = atoi(argv[++i]);
		i = (i == argc-1) ? (i) : (i+1);
		switch (model_number) {

			case 1: {
				if (string(argv[i]) == "--arg") {
					vector<float> model_param = cstring_to_string (i+1, argc-1, argv);
					float rho = model_param[0];
					float phi = Math::deg2rad(model_param[1]);
					cout << "Athimus angle: " << phi << endl;
					ReadyModel::linear_Ex_from_ct(rho, phi, 0.5);
				}
				else ReadyModel::linear_Ex_from_ct(0, 0, 0.01);
				break;
			} 

			case 2: {
				if (string(argv[i]) == "--arg") {
					vector<float> model_param = cstring_to_string (i+1, argc-1, argv);
					float z = model_param[0];
					cout << "Distance from source: " << z << endl;
					ReadyModel::linear_Hy_from_ct(0,0,z);
				}
				else ReadyModel::linear_Hy_from_ct(0,0,0.01);
				break;
			} 

			case 3: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::linear_Hy_from_ct_z();
				break;
			} 

			case 4: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::linear_Ex_from_ct_z();
				break;
			} 

			case 5: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::linear_Ex_from_rho_phi();
				break;
			} 

			case 6: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::linear_Ex_from_ct_rho();
				break;
			} 

			case 7: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::missile_effect_length();
				break;
			} 

			case 8: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::magnetic_static_magnitude();
				break;
			} 

			case 9: {
				if (string(argv[i]) == "--arg") {
					cerr << "--arg not allowed here" << endl;
					return -1;
				}
				ReadyModel::magnetic_static_terms();
				break;
			} 

			case 10: {
				ReadyModel::inegration_compare_i1();
				break;
			} 

			case 11: {
				ReadyModel::inegration_compare_i2();
				break;
			} 

			default:
				cerr << "Illegal model number: ";
				cerr << argv[i] << endl;
				return -1;
		}
	}

	/* else {
		cerr << "Illegal argument: ";
		cerr << argv[i] << endl;	
		cerr << "See --help for more information" << endl;			
		return -1;
	} */

	return 0;
}

vector<float> cstring_to_string (size_t from, size_t to, char* argv[])
{
	vector<float> res;
	for (size_t iter = from; iter <= to; iter++) {
		float val = atof(argv[iter]);
		res.push_back(val);
	}
	return res;
}

void print_configuration ()
{
	cout << boolalpha;
	cout << "Current configuration state ========================" << endl;
	cout << "Display progress: " << Config::print_progress() << endl;
	cout << "Selected color map: " << Config::plot_color_map() << endl;
	cout << "Display plot grid: " << Config::plot_grid() << endl;
	cout << "Display plot bound cadge: " << Config::plot_baund_cage() << endl;
	cout << "Selected plot device: " << Config::plot_device() << endl;
	cout << "Gnuplot binary: " << Config::path_gnuplot_binary() << endl;
	cout << "GNP execution script: " << Config::gnp_script_path() << endl;
	cout << "Magnetic field term number: " << Config::magnetic_term_num() << endl;
	cout << "Float bitrate for GMP: " << Config::float_bitrate() << endl;
	cout << "====================================================" << endl;
}
