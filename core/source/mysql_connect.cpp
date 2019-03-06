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
const std::string MySQL::SELECT_RESULT     = "SELECT result FROM maxwell.maxwell_data WHERE id = $POINT;";

MySQL::MySQL (const std::string& host, const std::string& user, const std::string& pass, const std::string& db)
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
}

MySQL::~MySQL ()
{
	mysql_close(this->connection);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	delete this->connection;
}

std::vector<std::pair<std::size_t,std::string>> MySQL::get_saved_problem_list () const
{
	int error_code = mysql_query(this->connection, MySQL::LIST_SAVED_MODELS.c_str());
	MySQL::throw_error_code(error_code);
	MYSQL_RES* serch_result = mysql_store_result(connection);
	std::size_t problems_num = mysql_num_rows(serch_result);

	std::vector<std::pair<std::size_t,std::string>> res;
	for (std::size_t item = 0; item < problems_num; item++) {
		MYSQL_ROW serch_row = mysql_fetch_row(serch_result);
		res.emplace_back(std::stoi(serch_row[0]),serch_row[1]);
	}

	mysql_free_result(serch_result);
	return res;
}

std::size_t MySQL::get_selected_problem () const
{
	if (!this->problem_selected) throw std::logic_error("Problem is not selected!");
	return this->problem_id;
}

bool MySQL::select_problem (std::size_t id)
{

	std::string problem_exists = MySQL::PRIBLEM_EXISTS;
	problem_exists = std::regex_replace(problem_exists, std::regex("\\$PROBLEM"), std::to_string(id));

	auto error_code = mysql_query(this->connection, problem_exists.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* result_exist = mysql_store_result(connection);
	bool exists = std::stoi(mysql_fetch_row(result_exist)[0]);
	mysql_free_result(result_exist);

	if (exists) {
		this->problem_selected = true;
		this->problem_id = id;
		return true;
	}

	return false;
}

bool MySQL::insert_problem (std::size_t id, const std::string& comment)
{
	std::string request = MySQL::INSERT_PROBLEM;
	request = std::regex_replace(request, std::regex("\\$COMMENT"), comment);
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(id));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	return true;
}

bool MySQL::select_point (double ct, double rho, double phi, double z)
{
	if (!problem_selected) throw std::logic_error("No problem_id selected");
	this->point_selected = true;

	std::string request = MySQL::SELECT_POINT;
	request = std::regex_replace(request, std::regex("\\$HEAD"), std::to_string(problem_id));
	request = std::regex_replace(request, std::regex("\\$TIME"), std::to_string(ct));
	request = std::regex_replace(request, std::regex("\\$RADIAL"), std::to_string(rho));
	request = std::regex_replace(request, std::regex("\\$AZIMUTH"), std::to_string(phi));
	request = std::regex_replace(request, std::regex("\\$DISTANCE"), std::to_string(z));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(this->connection);
	MYSQL_ROW row = mysql_fetch_row(responce);

	if (row != NULL) {
		this->point_selected = true;
		this->point_id = std::stod(row[0]);
		this->result = row[1] ? std::stod(row[1]) : NAN;
		mysql_free_result(responce);
		return true;
	}

	return false;
}

bool MySQL::insert_point (double ct, double rho, double phi, double z)
{
	if (!problem_selected) return false;

	std::string request = MySQL::INSERT_POINT;
	request = std::regex_replace(request, std::regex("\\$HEAD"), std::to_string(problem_id));
	request = std::regex_replace(request, std::regex("\\$TIME"), std::to_string(ct));
	request = std::regex_replace(request, std::regex("\\$RADIAL"), std::to_string(rho));
	request = std::regex_replace(request, std::regex("\\$AZIMUTH"), std::to_string(phi));
	request = std::regex_replace(request, std::regex("\\$DISTANCE"), std::to_string(z));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	this->point_id = mysql_insert_id(this->connection);
	this->result = NAN;
	this->point_selected = true;
	return true;
}

std::string MySQL::get_hostname() const
{
	return this->username + "@" + this->hostname;
}

bool MySQL::reconnect () const
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

double MySQL::get_result () const
{
	if (!point_selected || !problem_selected)
		throw std::logic_error("problem_id or point_id is not selected");
	return this->result;
}

bool MySQL::update_result (double value)
{
	std::string update_result = MySQL::UPDATE_RESULT;
	update_result = std::regex_replace(update_result, std::regex("\\$POINT"), std::to_string(point_id));
	update_result = std::regex_replace(update_result, std::regex("\\$VALUE"), std::to_string(value));

	int error_code = mysql_query(this->connection, update_result.c_str());
	MySQL::throw_error_code(error_code);
	this->result = value;
	return true;
}
