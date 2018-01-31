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

struct MySQL {

	MySQL (Config* global_config);
	std::string get_hostname() const;
	void reset_noise () const;
	void select_point (double ct, double rho, double phi, double z);
	~MySQL ();

	double get_noise () const;
	void set_noise (double value);

	double get_linear () const;
	void set_linear (double value);

	double get_square () const;
	void set_square (double value);

	double get_kerr () const;
	void set_kerr (double value);

	double get_value(const std::type_info&) const;
	void set_value(const std::type_info&, double value);

protected:

	static void throw_error_code (int code);
	static std::string to_string(const FieldComponent& type);
	void reconnect () const;

private:

	std::size_t problem_id;
	std::size_t point_id;
	Config* global_config;
	MYSQL* connection;

	double noise;
	double linear;
	double square;
	double kerr;

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
};

#endif /* mysql_connect_hpp */
