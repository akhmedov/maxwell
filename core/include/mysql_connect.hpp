//
//  mysql_connect.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 24.01.18.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef mysql_connect_hpp
#define mysql_connect_hpp

#include <mysql.h>

#include <thread>
#include <chrono>

#include <typeinfo> // typeid()
#include <cmath> // NAN
#include <regex> // std::replace()
#include <vector>
#include <string>
#include <exception>

struct problem_data { 
	std::size_t problem; // maxwell.maxwell_header.ID
	std::size_t point; // maxwell.maxwell_data.ID
	double result; // result
};

struct MySQL {

	MySQL (std::string host, std::string user, std::string pass, std::string db, std::vector<int> problem_id);
	~MySQL ();

	std::string get_hostname() const;
	void select_point (double ct, double rho, double phi, double z);
	double get_result (std::size_t problem) const;
	void set_result (std::size_t problem, double value);

protected:

	static void throw_error_code (int code);
	void reconnect () const;

private:

	MYSQL* connection;

	std::string hostname;
	std::string username;
	std::string password;
	std::string database;

	std::vector<problem_data> working;

	static const std::string USE_MAXWELL;
	static const std::string SET_WAIT_TIMEOUT;

	static const std::string DELETE_PROBLEM;
	static const std::string LIST_SAVED_MODELS;
	static const std::string INSERT_PROBLEM;
	static const std::string UPDATE_COMMENT;
	static const std::string PRIBLEM_EXISTS;

	static const std::string SELECT_POINT;
	static const std::string INSERT_POINT;
	static const std::string UPDATE_RESULT;
};

#endif /* mysql_connect_hpp */
