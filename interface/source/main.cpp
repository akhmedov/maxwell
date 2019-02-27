//
//  main.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include <string> // compare substr
#include <vector> // push_back

#include "cl_interface.hpp"
#include "plot_model.hpp"
#include "config.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	Config* global_conf = new Config();
	CLI* interface = new CLI(global_conf);

	PlotModel* plot_model = new PlotModel(global_conf);
	interface->set_mvc_model(plot_model);

	/* DataModel* data_model = new DataModel(global_conf);
	interface->set_mvc_model(data_model); */

	vector<string> cli_options(argv, argv + argc);
	interface->call_handler(cli_options);

	delete plot_model;
	delete interface;
	return 0;
}