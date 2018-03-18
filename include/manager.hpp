//
//  manager.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 23.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef manager_hpp
#define manager_hpp

#include "mysql_connect.hpp"
#include "abstract_field.hpp"
#include "uniform_disk_current.hpp"

#include <typeinfo>
#include <typeindex>

#include <set>
#include <mutex>
#include <stack>
#include <thread>
#include <chrono>
#include <functional>
#include <iostream>

typedef std::function<double(AbstractField*,double,double,double,double)> Component;

using std::chrono_literals::operator""s;
using std::chrono_literals::operator""ms;

struct TimeSort {
	bool operator() (const std::vector<double>& l_wp, const std::vector<double>& r_wp) const 
	{ return l_wp[0] < r_wp[0]; }
};

struct Manager {
	
	Manager ();
	Manager (std::size_t threads);
	void progress_bar (bool status = true);
	void add_argument (std::vector<double> argument);
	virtual std::vector<std::vector<double>> get_value ();

	void call ( std::function<double(double)> );
	void call ( std::vector<std::function<double(double)>> funcs);
	void call ( std::function<double(double,double,double,double)> );
	void call ( std::function<double(AbstractField*,double,double,double,double)>, AbstractField*);
	virtual void call ( std::vector<std::pair<Component,AbstractField*>>);

protected:
	void reset ();
	bool is_ready ();

	std::vector< std::thread > thread_list;
	std::vector< std::stack< std::vector<double> > > argument;
	std::multiset< std::vector<double>, TimeSort > result;
	std::size_t thread_number;
	std::mutex result_write;
	std::size_t data_left;
	std::size_t total_data;
	bool print_progtess = false;
};

struct SafeManager : public Manager {

	// SafeManger (Config* global);
	SafeManager (std::size_t threads, Config* global);
	std::vector<std::vector<double>> get_value ();
	void call ( std::vector<std::pair<Component,AbstractField*>>);
	~SafeManager ();

private:

	std::vector<MySQL*> client;
};

#endif /* manager_hpp */
