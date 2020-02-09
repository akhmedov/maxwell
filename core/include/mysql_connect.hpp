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

#include <utility>
#include <cmath> // NAN
#include <regex> // std::replace()
#include <vector>
#include <string>
#include <exception>

struct MySQL {

	MySQL (const std::string& host, const std::string& user, const std::string& pass, const std::string& db);
	~MySQL ();

	std::string get_hostname() const;
	std::vector<std::pair<std::size_t,std::string>> get_saved_problems () const;

	void delete_problem ();
	void select_problem (std::size_t id);
	void update_problem (const std::string& comment);
	std::size_t insert_problem (const std::string& comment);

	void delete_point ();
	void select_point (double ct, double rho, double phi, double z);
	std::size_t insert_point (double ct, double rho, double phi, double z);
	std::vector<std::size_t> select_entity ();
	double get_result () const;
	void update_result (double value);

protected:

	static void throw_error_code (int code);
	void reconnect () const;

private:

	MYSQL* connection;

	std::string hostname;
	std::string username;
	std::string password;
	std::string database;

	std::size_t problem_id {}; // maxwell.maxwell_header.ID
	std::size_t point_id {};   // maxwell.maxwell_data.ID
	double result {NAN};

	bool problem_selected {};
	bool point_selected {};

	static const std::string USE_MAXWELL;
	static const std::string SET_WAIT_TIMEOUT;

	static const std::string DELETE_PROBLEM;
	static const std::string LIST_SAVED_MODELS;
	static const std::string INSERT_PROBLEM;
	static const std::string UPDATE_COMMENT;
	static const std::string PRIBLEM_EXISTS;
	static const std::string SELECT_ENTITY;

	static const std::string DELETE_OBSERVER;
	static const std::string SELECT_OBSERVER;
	static const std::string INSERT_OBSERVER;
	static const std::string UPDATE_RESULT;
	static const std::string SELECT_RESULT;

	static const std::string DELETE_MODE;
	static const std::string SELECT_MODE;
	static const std::string INSERT_MODE;
	static const std::string UPDATE_EVOLUTION;
	static const std::string SELECT_EVOLUTION;
};

#endif /* mysql_connect_hpp */
