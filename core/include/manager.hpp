//
//  manager.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 23.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef manager_hpp
#define manager_hpp

#include "logger.hpp"
#include "mysql_connect.hpp"
#include "abstract_field.hpp"å

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
typedef std::function<double(AbstractField*,double,double,double,double,double)> Energy;

using std::chrono_literals::operator""s;
using std::chrono_literals::operator""ms;

template <std::size_t N> struct Sort {
	bool operator() (const std::vector<double>& l_wp, const std::vector<double>& r_wp) const 
	{ return l_wp[N] < r_wp[N]; }
};

template <std::size_t N> struct Manager {
	
	Manager (Logger* global_logger = NULL);
	Manager (std::size_t threads, Logger* global_logger = NULL);
	void progress_bar (bool status = true);
	void add_argument (std::vector<double> argument);
	virtual std::vector<std::vector<double>> get_value ();

	void call ( std::function<double(double)> );
	virtual void call ( std::vector<std::pair<Component,AbstractField*>>);
	virtual void call ( std::vector<std::pair<Energy,AbstractField*>>);

protected:
	void reset ();
	bool is_ready ();

	std::mutex active_thread;
	std::vector<bool> is_active;
	std::vector< std::thread > thread_list;
	std::vector< std::stack< std::vector<double> > > argument;
	std::multiset< std::vector<double>, Sort<N> > result;
	std::size_t thread_number;
	std::mutex result_write;
	std::size_t data_left;
	std::size_t total_data;
	bool print_progtess = false;

private:
	Logger* global_logger;
};

template <std::size_t N> struct SafeManager : public Manager<N> {

	// SafeManger (Config* global, Logger* global_logger = NULL);
	SafeManager (std::size_t threads, Config* global, Logger* global_logger = NULL);
	SafeManager (std::size_t threads, std::vector<Config*> gl_config, Logger* logger_ptr = NULL);
	std::vector<std::vector<double>> get_value ();
	void call ( std::vector<std::pair<Component,AbstractField*>>);
	void call ( std::vector<std::pair<Energy,AbstractField*>>);
	~SafeManager ();

private:
	const std::size_t problems;
	std::vector<MySQL*> client;
};

//=============================================================================
//== Manger ===================================================================
//=============================================================================

template <std::size_t N> Manager<N>::Manager (Logger* logger_ptr)
: Manager(std::thread::hardware_concurrency() - 1, logger_ptr) 
{ }

//=============================================================================

template <std::size_t N> Manager<N>::Manager (std::size_t threads, Logger* logger_ptr)
{
	this->global_logger = logger_ptr;
	this->data_left = 0;
	this->total_data = 0;
	this->thread_number = threads;
	for (std::size_t iter = 0; iter < this->thread_number; iter++) {
		this->argument.push_back( std::stack<std::vector<double>>() );
		this->is_active.push_back(false);
	}
}

//=============================================================================

template <std::size_t N> void Manager<N>::reset ()
{
	// TODO: kill chiled threads
	this->data_left = 0;
	this->total_data = 0;
	this->argument.clear();
	this->result.clear();
}

//=============================================================================

template <std::size_t N> void Manager<N>::progress_bar (bool status)
{
	this->print_progtess = status;
}

//=============================================================================

template <std::size_t N> std::vector<std::vector<double>> Manager<N>::get_value ()
{
	while (!this->is_ready()) { }

	auto tmp = this->result;
	std::vector<std::vector<double>> res(tmp.size());
	std::copy(tmp.begin(), tmp.end(), res.begin());

	if (this->global_logger) this->global_logger->write_all();
	
	this->reset();
	return res;
}

//=============================================================================


template <std::size_t N> bool Manager<N>::is_ready ()
{
	bool ready = true;

	for (auto&& i : this->argument) 
		if (!i.empty()) ready = false;
	for (auto&& i : this->thread_list) 
		if (i.joinable()) ready = false;
	for (auto&& i : this->is_active)
		if (i) ready = false;

	return ready;
}

//=============================================================================

template <std::size_t N> void Manager<N>::add_argument (std::vector<double> arg)
{
	this->data_left++;
	this->total_data++;
	static std::size_t thread = 0;
	this->argument[thread].push(arg);
	if (++thread == this->thread_number) thread = 0;
}

//=============================================================================

template <std::size_t N> void Manager<N>::call ( std::function<double(double)> func)
{
	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			this->active_thread.lock();
			this->is_active[i] = true;
			this->active_thread.unlock();
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();
			double res = func(arg[0]);
			arg.insert( arg.begin(), res );
			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
			this->active_thread.lock();
			this->is_active[i] = false;
			this->active_thread.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	std::size_t last_printed = 101;
	while (!this->is_ready()) {
		if (this->print_progtess) {
			double progress = (double) this->data_left / (double) this->total_data;
			progress *= 100;
			if ((std::size_t)progress != last_printed) {
				last_printed = (std::size_t) progress;
				std::cout << '\r' << "Evaluation progress: " << 100 - (std::size_t) progress << '%';
				std::cout.flush();
			}
		}
	}
	std::cout << std::endl;
}

//=============================================================================

template <std::size_t N> void Manager<N>::call ( std::vector<std::pair<Component,AbstractField*>> field )
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			this->active_thread.lock();
			this->is_active[i] = true;
			this->active_thread.unlock();
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();

			for (auto f : field) {
				double res = f.first(f.second, arg[0], arg[1], arg[2], arg[3]);
				if (std::isnan(res)) throw std::logic_error("Error: NAN model value.");
				arg.push_back(res);
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			this->result_write.unlock();
			this->active_thread.lock();
			this->is_active[i] = false;
			this->active_thread.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	std::size_t last_printed = 101;
	while (!this->is_ready()) {
		if (this->print_progtess) {
			double progress = (double) this->data_left / (double) this->total_data;
			progress *= 100;
			if ((std::size_t)progress != last_printed) {
				last_printed = (std::size_t) progress;
				std::cout << '\r' << "Evaluation progress: " << 100 - (std::size_t) progress << '%';
				std::cout.flush();
			}
		}
	}
	std::cout << std::endl;
}

//=============================================================================

template <std::size_t N> void Manager<N>::call ( std::vector<std::pair<Energy,AbstractField*>> field )
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			this->active_thread.lock();
			this->is_active[i] = true;
			this->active_thread.unlock();
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();

			for (auto f : field) {
				double res = f.first(f.second, arg[0], arg[1], arg[2], arg[3], arg[4]);
				if (std::isnan(res)) res = 0; // throw std::logic_error("Error: NAN model value.");
				arg.push_back(res);
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
			this->active_thread.lock();
			this->is_active[i] = false;
			this->active_thread.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	std::size_t last_printed = 101;
	while (!this->is_ready()) {
		if (this->print_progtess) {
			double progress = (double) this->data_left / (double) this->total_data;
			progress *= 100;
			if ((std::size_t)progress != last_printed) {
				last_printed = (std::size_t) progress;
				std::cout << '\r' << "Evaluation progress: " << 100 - (std::size_t) progress << '%';
				std::cout.flush();
			}
		}
	}
	std::cout << std::endl;
}

//============================================================================================================
//== SafeManger ==============================================================================================
//============================================================================================================

/* SafeManager::SafeManager (Config* gl_config, Logger* logger_ptr) 
: Manager(), client(gl_config) { } */

template <std::size_t N> SafeManager<N>::SafeManager (std::size_t threads, std::vector<Config*> gl_config, Logger* logger_ptr) 
: Manager<N>(threads,logger_ptr), problems(1), client()
{
	for (std::size_t c = 0; c < this->thread_number; c++) {
		MySQL* thread_client = new MySQL(gl_config);
		this->client.push_back(thread_client);
	}
}

template <std::size_t N> SafeManager<N>::SafeManager (std::size_t threads, Config* gl_config, Logger* logger_ptr) 
: Manager<N>(threads,logger_ptr), problems(1), client()
{
	for (std::size_t c = 0; c < this->thread_number; c++) {
		MySQL* thread_client = new MySQL(gl_config);
		this->client.push_back(thread_client);
	}
}

//============================================================================================================

template <std::size_t N> void SafeManager<N>::call (std::vector<std::pair<Component,AbstractField*>> field)
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			this->active_thread.lock();
			this->is_active[i] = true;
			this->active_thread.unlock();
			std::vector<double> arg = this->argument[i].top();
			this->argument[i].pop();
			
			this->client[i]->select_point(arg[0], arg[1], arg[2], arg[3]);

			for (std::size_t j = 0; j < this->problems; j++) {
				for (auto f : field) {

					const std::type_info& type = typeid(*f.second);
					double db_res = this->client[i]->get_value(j,type);
					double res;

					if ( std::isnan(db_res) ) {
						res = f.first(f.second, arg[0], arg[1], arg[2], arg[3]);
						if (std::isnan(res)) throw std::logic_error("Error: NAN model value.");
						this->client[i]->set_value(j,type,res);
					} else {
						res = db_res;
					}

					arg.push_back(res);
				}
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			this->result_write.unlock();
			this->active_thread.lock();
			this->is_active[i] = false;
			this->active_thread.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	std::size_t last_printed = 101;
	while (!this->is_ready()) {
		if (this->print_progtess) {
			double progress = (double) this->data_left / (double) this->total_data;
			progress *= 100;
			if ((std::size_t)progress != last_printed) {
				last_printed = (std::size_t) progress;
				std::cout << '\r' << "Evaluation progress: " << 100 - (std::size_t) progress << '%';
				std::cout.flush();
			}
		}
	}
	std::cout << std::endl;
}

//============================================================================================================

template <std::size_t N> void SafeManager<N>::call ( std::vector<std::pair<Energy,AbstractField*>> field )
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			this->active_thread.lock();
			this->is_active[i] = true;
			this->active_thread.unlock();
			std::vector<double> arg = this->argument[i].top();
			this->argument[i].pop();

			this->client[i]->select_point(0, arg[0], arg[1], arg[2]);

			for (std::size_t j = 0; j < this->problems; j++) {
				for (auto f : field) {

					const std::type_info& type = typeid(*f.second);
					double db_res = this->client[i]->get_value(j,type);
					double res;

					if ( std::isnan(db_res) ) {
						res = f.first(f.second, arg[0], arg[1], arg[2], arg[3], arg[4]);
						if (std::isnan(res)) throw std::logic_error("Error: NAN model value.");
						this->client[i]->set_value(j,type,res);
					} else {
						res = db_res;
					}

					arg.push_back(res);
				}
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			this->result_write.unlock();
			this->active_thread.lock();
			this->is_active[i] = false;
			this->active_thread.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	std::size_t last_printed = 101;
	while (!this->is_ready()) {
		if (this->print_progtess) {
			double progress = (double) this->data_left / (double) this->total_data;
			progress *= 100;
			if ((std::size_t)progress != last_printed) {
				last_printed = (std::size_t) progress;
				std::cout << '\r' << "Evaluation progress: " << 100 - (std::size_t) progress << '%';
				std::cout.flush();
			}
		}
	}
	std::cout << std::endl;
}

//============================================================================================================

template <std::size_t N> std::vector<std::vector<double>> SafeManager<N>::get_value ()
{
	while (!this->is_ready()) { }

	auto tmp = this->result;
	std::vector<std::vector<double>> res(tmp.size());
	std::copy(tmp.begin(), tmp.end(), res.begin());

	// next section is not optimum - use for_all<for_all>

	for (auto line : res)
		for (auto item : line)
			if (std::isnan(item)) item = 0;
	
	this->reset();
	return res;
}

//============================================================================================================

template <std::size_t N> SafeManager<N>::~SafeManager ()
{
	for (auto i : this->client) delete i;
	this->client.clear();
}

#endif /* manager_hpp */
