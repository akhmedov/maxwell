//
//  mysql_connect.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 24.01.18.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "mysql_connect.hpp"

const std::string MySQL::SET_WAIT_TIMEOUT  = "SET wait_timeout=9999999;";
const std::string MySQL::USE_MAXWELL       = "USE maxwell;";

const std::string MySQL::LIST_SAVED_MODELS = "SELECT * FROM maxwell.maxwell_header;";
const std::string MySQL::INSERT_PROBLEM    = "INSERT INTO maxwell.maxwell_header SET id = $PROBLEM_ID, comment = $COMMENT;";
const std::string MySQL::UPDATE_COMMENT    = "UPDATE maxwell.maxwell_header SET comment = $COMMENT WHERE id = $PROBLEM;";
const std::string MySQL::DELETE_PROBLEM    = "DELETE FROM maxwell_header WHERE id = $PROBLEM;";
const std::string MySQL::PRIBLEM_EXISTS    = "SELECT EXISTS(SELECT id FROM maxwell.maxwell_header WHERE id = $PROBLEM);";

const std::string MySQL::SELECT_POINT      = "SELECT id,lineary,square,kerr FROM maxwell.maxwell_data WHERE head_id = $HEAD AND ct = $TIME AND rho = $RADIAL AND phi = $AZIMUTH AND z = $DISTANCE;";
const std::string MySQL::INSERT_POINT      = "INSERT INTO maxwell.maxwell_data SET head_id = $HEAD, ct = $TIME, rho = $RADIAL, phi = $AZIMUTH, z = $DISTANCE;";
const std::string MySQL::UPDATE_RESULT     = "UPDATE maxwell.maxwell_data SET lineary = $VALUE WHERE id = $POINT;";

MySQL::MySQL (std::string host, std::string user, std::string pass, std::string db, std::vector<int> problem_id)
: hostname(host), username(user), password(pass), database(db) 
{
	this->connection = mysql_init(NULL);

	auto connected = mysql_real_connect(this->connection,
										this->hostname.c_str(),
										this->username.c_str(),
										this->password.c_str(),
										this->database.c_str(),
										0, NULL, 0);

	if (!connected) throw std::logic_error("Connection failed! MySQL is not connected.");

	int error_code = mysql_query(this->connection, MySQL::SET_WAIT_TIMEOUT.c_str());
	MySQL::throw_error_code(error_code);

	error_code = mysql_query(this->connection, MySQL::USE_MAXWELL.c_str());
	MySQL::throw_error_code(error_code);

	for (auto id : problem_id) {

		this->working.push_back({(std::size_t)id, 0, 0});

		std::string problem_exists = MySQL::PRIBLEM_EXISTS;
		problem_exists = std::regex_replace(problem_exists, std::regex("\\$PROBLEM"), std::to_string(problem_id[0]));

		error_code = mysql_query(this->connection, problem_exists.c_str());
		MySQL::throw_error_code(error_code);

		MYSQL_RES* result_exist = mysql_store_result(connection);
		bool exists = mysql_fetch_row(result_exist)[0];
		mysql_free_result(result_exist);

		if (!exists) {

			std::string insert_problem = MySQL::INSERT_PROBLEM;
			insert_problem = std::regex_replace(insert_problem, std::regex("\\$COMMENT"), "no comment");
			insert_problem = std::regex_replace(insert_problem, std::regex("\\$PROBLEM"), std::to_string(id));

			error_code = mysql_query(this->connection, insert_problem.c_str());
			MySQL::throw_error_code(error_code);
		}
	}
}

void MySQL::select_point (double ct, double rho, double phi, double z)
{
	for (auto item : working) {

		std::string select_point = MySQL::SELECT_POINT;
		select_point = std::regex_replace(select_point, std::regex("\\$HEAD"), std::to_string(item.problem));
		select_point = std::regex_replace(select_point, std::regex("\\$TIME"), std::to_string(ct));
		select_point = std::regex_replace(select_point, std::regex("\\$RADIAL"), std::to_string(rho));
		select_point = std::regex_replace(select_point, std::regex("\\$AZIMUTH"), std::to_string(phi));
		select_point = std::regex_replace(select_point, std::regex("\\$DISTANCE"), std::to_string(z));

		int error_code = mysql_query(this->connection, select_point.c_str());
		MySQL::throw_error_code(error_code);

		MYSQL_RES* serch_result = mysql_store_result(this->connection);
		MYSQL_ROW serch_row = mysql_fetch_row(serch_result);

		if (serch_row != NULL) {

			item.point = std::stod(serch_row[0]);
			item.result = serch_row[1] ? std::stod(serch_row[1]) : NAN;
			mysql_free_result(serch_result);

		} else {

			std::string insert_point = MySQL::INSERT_POINT;
			insert_point = std::regex_replace(insert_point, std::regex("\\$HEAD"), std::to_string(item.problem));
			insert_point = std::regex_replace(insert_point, std::regex("\\$TIME"), std::to_string(ct));
			insert_point = std::regex_replace(insert_point, std::regex("\\$RADIAL"), std::to_string(rho));
			insert_point = std::regex_replace(insert_point, std::regex("\\$AZIMUTH"), std::to_string(phi));
			insert_point = std::regex_replace(insert_point, std::regex("\\$DISTANCE"), std::to_string(z));

			error_code = mysql_query(this->connection, insert_point.c_str());
			MySQL::throw_error_code(error_code);
			mysql_free_result(serch_result);

			item.point = mysql_insert_id(this->connection);
			item.result = NAN;
		}
	}
}

std::string MySQL::get_hostname() const
{
	return this->username + "@" + this->hostname;
}

void MySQL::reconnect () const
{
	throw std::logic_error("MySQL::reconnect not implemented");
}

void MySQL::throw_error_code (int code)
{
	std::this_thread::sleep_for(std::chrono::milliseconds(20));
	std::string message = "Connection esteblished. Error code: ";
	switch (code) {
		case 0: return;
		case 1: message += "no respound data"; break;
		// TODO: more codes
		default: message += std::to_string(code);
	}
	throw std::logic_error(message);
}

MySQL::~MySQL ()
{
	mysql_close(this->connection);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

//=======================================================================================

double MySQL::get_result (std::size_t problem) const
{
	return this->working[problem].result;
}

void MySQL::set_result (std::size_t problem, double value)
{
	this->working[problem].result = value;
	std::size_t point_id = this->working[problem].point;
	
	std::string update_result = MySQL::UPDATE_RESULT;
	update_result = std::regex_replace(update_result, std::regex("\\$POINT"), std::to_string(point_id));
	update_result = std::regex_replace(update_result, std::regex("\\$VALUE"), std::to_string(value));
	
	int error_code = mysql_query(this->connection, update_result.c_str());
	MySQL::throw_error_code(error_code);
}
