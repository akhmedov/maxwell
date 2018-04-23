//
//  manager.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 23.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "manager.hpp"

Manager::Manager (Logger* logger_ptr)
: Manager(std::thread::hardware_concurrency() - 1, logger_ptr) 
{ }

//=============================================================================

Manager::Manager (std::size_t threads, Logger* logger_ptr)
{
	this->global_logger = logger_ptr;
	this->data_left = 0;
	this->total_data = 0;
	this->thread_number = threads;
	for (std::size_t iter = 0; iter < this->thread_number; iter++)
		this->argument.push_back( std::stack<std::vector<double>>() );
}

//=============================================================================

void Manager::reset ()
{
	// TODO: kill chiled threads
	this->data_left = 0;
	this->total_data = 0;
	this->argument.clear();
	this->result.clear();
}

//=============================================================================

void Manager::progress_bar (bool status)
{
	this->print_progtess = status;
}

//=============================================================================

std::vector<std::vector<double>> Manager::get_value ()
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

bool Manager::is_ready ()
{
	bool ready = true;

	for (auto&& i : this->argument) 
		if (!i.empty()) ready = false;
	for (auto&& i : this->thread_list) 
		if (i.joinable()) ready = false;

	return ready;
}

//=============================================================================

void Manager::add_argument (std::vector<double> arg)
{
	this->data_left++;
	this->total_data++;
	static std::size_t thread = 0;
	this->argument[thread].push(arg);
	if (++thread == this->thread_number) thread = 0;
}

//=============================================================================

void Manager::call ( std::function<double(double)> func)
{
	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();
			double res = func(arg[0]);
			arg.insert( arg.begin(), res );
			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
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

void Manager::call ( std::function<double(double,double,double,double)> func)
{

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {

			std::vector<double> arg = this->argument[i].top();
			this->argument[i].pop();

			double res = func(arg[0], arg[1], arg[2], arg[3]);
			arg.insert( arg.begin(), res );

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			this->result_write.unlock();
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

void Manager::call ( std::vector<std::function<double(double)>> funcs)
{

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();

			for (auto f : funcs) {
				double res = f(arg.front());
				arg.push_back(res);
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
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

void Manager::call ( std::function<double(AbstractField*,double,double,double,double)> component, AbstractField* field)
{
	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();
			double res = component(field, arg[0], arg[1], arg[2], arg[3]);
			if (std::isnan(res)) throw std::logic_error("Error: NAN model value.");
			arg.push_back(res);
			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
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
	// std::this_thread::sleep_for(3s);
}

//=============================================================================

void Manager::call ( std::vector<std::pair<Component,AbstractField*>> field )
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
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
			result_write.unlock();
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
//== SafeManger ===============================================================
//=============================================================================

/* SafeManager::SafeManager (Config* gl_config, Logger* logger_ptr) 
: Manager(), client(gl_config) { } */

SafeManager::SafeManager (std::size_t threads, Config* gl_config, Logger* logger_ptr) 
: Manager(threads,logger_ptr), client()
{
	for (std::size_t c = 0; c < this->thread_number; c++) {
		MySQL* thread_client = new MySQL(gl_config);
		this->client.push_back(thread_client);
	}
}

void SafeManager::call (std::vector<std::pair<Component,AbstractField*>> field)
{
	if (!field.size()) throw std::invalid_argument("Manager::call argument must not be empty");

	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();
			
			this->client[i]->select_point(arg[0], arg[1], arg[2], arg[3]);

			for (auto f : field) {

				const std::type_info& type = typeid(*f.second);
				double db_res = this->client[i]->get_value(type);
				double res;

				if ( std::isnan(db_res) ) {
					res = f.first(f.second, arg[0], arg[1], arg[2], arg[3]);
					if (std::isnan(res)) throw std::logic_error("Error: NAN model value.");
					this->client[i]->set_value(type,res);
				} else {
					res = db_res;
				}

				arg.push_back(res);
			}

			this->result_write.lock();
			this->data_left--;
			this->result.insert( arg );
			result_write.unlock();
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

std::vector<std::vector<double>> SafeManager::get_value ()
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

SafeManager::~SafeManager ()
{
	for (auto i : this->client) delete i;
	this->client.clear();
}
