//
//  manager.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 23.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef manager_hpp
#define manager_hpp

#include "abstract_field.hpp"
#include "uniform_disk_current.hpp"

#include <mutex>
#include <stack>
#include <thread>
#include <chrono>
#include <functional>
#include <iostream>

using std::chrono_literals::operator""s;
using std::chrono_literals::operator""ms;

struct Manager {
	Manager ();
	Manager (std::size_t threads);
	void progress_bar (bool status = true);
	void call ( double (*func) (double) );
	// void call ( std::vector<double (*func) (double)> );
	void call ( std::vector<std::function<double(double)>> funcs);
	void call ( double (*field) (double, double, double, double) );
	void call ( const AbstractField& field, std::string method );
	void add_argument (std::vector<double> argument);
	bool is_ready ();
	std::stack<std::vector<double>> get_value ();
	std::vector< std::thread > thread_list;
private:
	void call ( const AbstractField& field, std::size_t method );
	std::vector< std::stack< std::vector<double> > > argument;
	std::stack< std::vector<double> > result;
	std::size_t thread_number;
	std::mutex result_write;
	std::size_t data_left;
	std::size_t total_data;
	bool print_progtess = false;
};

#endif /* manager_hpp */
