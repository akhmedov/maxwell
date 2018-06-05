//
//  mysql_connect.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 24.01.18.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef mysql_connect_hpp
#define mysql_connect_hpp

#include "config.hpp"
#include "noise.hpp"
#include "updisk_meandr.hpp"
#include "linear_duhamel.hpp"
#include "uniform_disk_current.hpp"
#include "kerr_amendment.hpp"

#include <iostream>

#include <mysql.h>
#include <thread>
#include <chrono>
#include <typeinfo> // typeid()
#include <cmath> // NAN
#include <regex> // std::replace()
#include <vector>
#include <array>
#include <string>
#include <exception>

typedef std::pair<std::size_t,std::size_t> id_id;

struct MySQL {

	MySQL (Config* global_config);
	MySQL (const std::vector<Config*> &vect);
	
	std::string get_hostname() const;
	void select_point (double ct, double rho, double phi, double z);
	~MySQL ();

	double get_linear (std::size_t problem) const;
	void set_linear (std::size_t problem, double value);

	double get_square (std::size_t problem) const;
	void set_square (std::size_t problem, double value);

	double get_kerr (std::size_t problem) const;
	void set_kerr (std::size_t problem, double value);

	double get_value(std::size_t problem, const std::type_info&) const;
	void set_value(std::size_t problem, const std::type_info&, double value);

protected:

	static void throw_error_code (int code);
	static std::string to_string(const ImpulseShape& type);
	static std::string to_string(const FieldComponent& type);
	void reconnect () const;

private:

	std::vector<id_id> item;
	std::vector<Config*> global_conf;
	MYSQL* connection;

	// double noise;
	std::vector<double> linear;
	std::vector<double> square;
	std::vector<double> kerr;

	static const std::string USE_MAXWELL;
	static const std::string SELECT_PROBLEM_ID;
	static const std::string INSERT_PROBLEM;
	static const std::string SELECT_POINT;
	static const std::string INSERT_POINT;
	static const std::string UPDATE_NOISE;
	static const std::string UPDATE_LINEAR;
	static const std::string UPDATE_SQUARE;
	static const std::string UPDATE_KERR;
	static const std::string RESET_NOISE;
	static const std::string SET_WAIT_TIMEOUT;
};

#endif /* mysql_connect_hpp */
