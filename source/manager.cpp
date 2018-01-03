//
//  manager.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 23.07.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "manager.hpp"

Manager::Manager ()
{
	this->data_left = 0;
	this->total_data = 0;
	this->thread_number = std::thread::hardware_concurrency() - 1;
	for (std::size_t iter = 0; iter < this->thread_number; iter++)
		this->argument.push_back( std::stack<std::vector<double>>() );
}

//=============================================================================

Manager::Manager (std::size_t threads)
{
	this->thread_number = threads;
	for (std::size_t iter = 0; iter < this->thread_number; iter++)
		this->argument.push_back( std::stack<std::vector<double>>() );
}

//=============================================================================

void Manager::reset ()
{
	this->data_left = 0;
	this->total_data = 0;
	this->argument.clear();
	// this->result.clear();
}

//=============================================================================

void Manager::progress_bar (bool status)
{
	this->print_progtess = status;
}

//=============================================================================

std::stack<std::vector<double>> Manager::get_value ()
{
	return this->result;
}

//=============================================================================

bool Manager::is_ready ()
{
	bool ready = true;
	for (auto&& i : this->argument) ready = (ready && i.empty());
	return ready && this->result_write.try_lock();
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
			this->result.push( arg );
			result_write.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	// called = true;
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
			argument[i].pop();
			double res = func(arg[0], arg[1], arg[2], arg[3]);
			arg.insert( arg.begin(), res );
			this->result_write.lock();
			this->data_left--;
			this->result.push( arg );
			result_write.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	// called = true;
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
				// arg.insert( arg.begin(), res );
				arg.push_back(res);
			}

			this->result_write.lock();
			this->data_left--;
			this->result.push( arg );
			result_write.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	// called = true;
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

void Manager::call ( std::function<double(const AbstractField&,double,double,double,double)> component, const AbstractField& field)
{
	auto call_thread = [&] (std::size_t i) {
		while (!this->argument[i].empty()) {
			std::vector<double> arg = this->argument[i].top();
			argument[i].pop();
			double res = component(field, arg[0], arg[1], arg[2], arg[3]);
			arg.insert( arg.begin(), res );
			this->result_write.lock();
			this->data_left--;
			this->result.push( arg );
			result_write.unlock();
		}
	};

	this->thread_list = std::vector<std::thread>(this->thread_number);

	std::size_t num = 0;
	for (auto& i : this->thread_list) {
		i = std::thread(call_thread, num++);
		i.detach();
	}

	// called = true;

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
