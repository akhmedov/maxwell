//
//  cl_interface.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 28.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef cl_interface_hpp
#define cl_interface_hpp

#include "noise.hpp"
#include "updisk_meandr.hpp"
#include "kerr_amendment.hpp"
#include "linear_duhamel.hpp"
#include "uniform_disk_current.hpp"

#include "config.hpp"
#include "plot_model.hpp"

#include <iostream>
#include <fstream> // cheack if file open, config file read
#include <string>
#include <sstream>
#include <vector>
#include <utility> // std::make_pair
#include <set>
#include <map>

struct CLI {
	CLI (Config* global_conf);
	void call_handler (const std::vector<std::string>& argv);
	void set_mvc_model(PlotModel* model);
	std::vector<std::pair<Component,AbstractField*>> em_problem (const ImpulseShape& problem_name, const Component& comp) const;
	~CLI ();

protected:
	void read_config_file () const;
	void update_config (const std::string& param, const std::string& arg) const;
	void update_config (const std::string& param) const;
	void run_list ();

	void print_version () const;
	void print_help () const;
	void print_arguments() const;

	static bool is_float (const std::string& literal);
	static bool is_int (const std::string& literal);
	static bool is_path (const std::string& literal);
	static bool is_pair (const std::string& literal);

	static bool is_onedim_model (const std::string& literal);
	static bool is_twodim_model (const std::string& literal);
	static bool is_const_value (const std::string& literal);
	static bool is_range_value (const std::string& literal);
	static bool is_event (const std::string& literal);
	
	static std::vector<std::string> split (const std::string& literal, char separator);
	static Component func_ptr_of (const FieldComponent& comp);

private:
	Logger* global_log;
	Config* global_conf;
	PlotModel* plot_model;
	std::set<std::string> unar_cmd;
	std::vector<std::string> argtype_regs;
	std::multimap<std::string,std::string> binar_cmd;
	std::multimap<std::string,std::string> binar_def_cmd;
	std::multimap<std::string,std::string> binar_var_cmd;
};

#endif /* cl_interface_hpp */
