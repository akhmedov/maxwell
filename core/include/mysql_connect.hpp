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

	void delete_probe ();
	void select_probe (double ct, double rho, double phi, double z);
	std::size_t insert_probe (double ct, double rho, double phi, double z);
	std::vector<std::size_t> select_all_probes (); // get all probe->id for selected problem->id
	double get_probe_result () const;
	void update_probe_result (double value);

	void delete_coeff ();
	void select_coeff (int m, double nu);
	std::size_t insert_coeff (int m, double nu);
	std::vector<std::size_t> select_all_coeffs (); // get all coeff->id for selected probe->id
	double get_coeff_result () const;
	void update_coeff_result (double value);

protected:

	static void throw_error_code (int code);
	void reconnect () const;

private:

	MYSQL* connection;

	std::string hostname;
	std::string username;
	std::string password;
	std::string database;

	std::size_t problem_id {}; // maxwell.problem.id
	std::size_t probe_id {};   // maxwell.probe.id
	std::size_t evo_id {};   // maxwell.evolution.id
	double result {NAN}; // maxwell.probe.result
	double vmh {NAN}; // maxwell.evolution.vmh

	bool problem_selected {};
	bool probe_selected {};
	bool evolution_selected {};

	static const std::string USE_MAXWELL;
	static const std::string SET_WAIT_TIMEOUT;

	static const std::string DELETE_PROBLEM;
	static const std::string LIST_SAVED_MODELS;
	static const std::string INSERT_PROBLEM;
	static const std::string UPDATE_COMMENT;
	static const std::string PRIBLEM_EXISTS;
	static const std::string SELECT_ALL_PROBES;

	static const std::string DELETE_PROBE;
	static const std::string SELECT_PROBE;
	static const std::string INSERT_PROBE;
	static const std::string UPDATE_PROBE_RESULT;
	static const std::string SELECT_PROBE_RESULT;

	static const std::string DELETE_COEFF;
	static const std::string SELECT_COEFF;
	static const std::string INSERT_COEFF;
	static const std::string UPDATE_COEFF_VAL;
	static const std::string SELECT_COEFF_VAL;
	static const std::string SELECT_ALL_COEFF;
};

#endif /* mysql_connect_hpp */
