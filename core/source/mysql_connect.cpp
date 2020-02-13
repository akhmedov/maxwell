//
//  mysql_connect.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 24.01.18.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "mysql_connect.hpp"

const std::string MySQL::SET_WAIT_TIMEOUT  		= "SET wait_timeout=9999999;";
const std::string MySQL::USE_MAXWELL       		= "USE maxwell;";

const std::string MySQL::LIST_SAVED_MODELS 		= "SELECT * FROM problem ORDER BY id;";
const std::string MySQL::INSERT_PROBLEM    		= "INSERT INTO problem SET comment = '$COMMENT';";
const std::string MySQL::UPDATE_COMMENT    		= "UPDATE problem SET comment = '$COMMENT' WHERE id = $PROBLEM;";
const std::string MySQL::DELETE_PROBLEM    		= "DELETE FROM problem WHERE id = $PROBLEM;";
const std::string MySQL::SELECT_PROBLEM_ID    	= "SELECT EXISTS(SELECT id FROM problem WHERE id = $PROBLEM);";
const std::string MySQL::SELECT_PROBLEM_COMMENT = "SELECT id FROM problem WHERE comment = '$COMMENT';";
const std::string MySQL::SELECT_ALL_PROBES		= "SELECT id FROM probe WHERE problem_id = $PROBLEM;";

const std::string MySQL::DELETE_PROBE   		= "DELETE FROM probe WHERE id = $POINT;";
const std::string MySQL::SELECT_PROBE   		= "SELECT id,result FROM probe WHERE problem_id = $HEAD AND ct = $TIME AND rho = $RADIAL AND phi = $AZIMUTH AND z = $DISTANCE;";
const std::string MySQL::INSERT_PROBE   		= "INSERT INTO probe SET problem_id = $HEAD, ct = $TIME, rho = $RADIAL, phi = $AZIMUTH, z = $DISTANCE;";
const std::string MySQL::UPDATE_PROBE_RESULT    = "UPDATE probe SET result = $VALUE WHERE id = $POINT;";
const std::string MySQL::SELECT_PROBE_RESULT    = "SELECT result FROM probe WHERE id = $POINT;";

const std::string MySQL::DELETE_COEFF	   		= "DELETE FROM evolution WHERE id = $INDEX;";
const std::string MySQL::SELECT_COEFF 	   		= "SELECT id,vmh FROM evolution WHERE problem_id = $PROBLEM AND probe_id = $PROBE AND m = $M AND nu = $NU;";
const std::string MySQL::INSERT_COEFF 	 		= "INSERT INTO evolution SET problem_id = $PROBLEM, probe_id = $PROBE, m = $M, nu = $NU;";
const std::string MySQL::UPDATE_COEFF_VAL  		= "UPDATE evolution SET vmh = $VALUE WHERE id = $INDEX;";
const std::string MySQL::SELECT_COEFF_VAL 		= "SELECT vmh FROM evolution WHERE id = $INDEX;";
const std::string MySQL::SELECT_ALL_COEFF  		= "SELECT nu,vmh FROM evolution WHERE problem_id = $PROBLEM AND probe_id = $PROBE $M_OPT AND vmh IS NOT NULL ORDER BY nu;";

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

	if (!connected) {
		// delete this->connection;
		throw std::logic_error("Connection failed! MySQL is not connected.");
	}

	int error_code = mysql_query(this->connection, MySQL::SET_WAIT_TIMEOUT.c_str());
	MySQL::throw_error_code(error_code);

	error_code = mysql_query(this->connection, MySQL::USE_MAXWELL.c_str());
	MySQL::throw_error_code(error_code);
}

MySQL::~MySQL ()
{
	mysql_close(this->connection);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	// delete this->connection;
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

std::vector<std::pair<std::size_t,std::string>> MySQL::get_saved_problems () const
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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void MySQL::select_problem (std::size_t id)
{
	std::string problem_exists = MySQL::SELECT_PROBLEM_ID;
	problem_exists = std::regex_replace(problem_exists, std::regex("\\$PROBLEM"), std::to_string(id));

	auto error_code = mysql_query(this->connection, problem_exists.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* result_exist = mysql_store_result(connection);
	bool exists = std::stoi(mysql_fetch_row(result_exist)[0]);
	mysql_free_result(result_exist);

	if (exists) this->problem_id = id;
	else this->problem_id = this->insert_problem("autoinserted problem [TIMESTAMP]");
	this->problem_selected = true;
}

void MySQL::select_problem (const std::string& comment)
{
	std::string request = MySQL::SELECT_PROBLEM_COMMENT;
	request = std::regex_replace(request, std::regex("\\$COMMENT"), comment);

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(this->connection);
	MYSQL_ROW row = mysql_fetch_row(responce);

	if (row != NULL) {
		this->problem_id = std::stod(row[0]);
		mysql_free_result(responce);
	} else {
		this->insert_problem(comment);
	}

	this->problem_selected = true;
}

void MySQL::delete_problem ()
{
	std::string request = MySQL::DELETE_PROBLEM;
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));
	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	this->probe_selected = false;
}

std::size_t MySQL::insert_problem (const std::string& comment)
{
	std::string request = MySQL::INSERT_PROBLEM;
	request = std::regex_replace(request, std::regex("\\$COMMENT"), comment);

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	this->problem_id = mysql_insert_id(this->connection);
	this->problem_selected = true;
	return mysql_insert_id(this->connection);
}

void MySQL::update_problem (const std::string& comment)
{
	std::string request = MySQL::UPDATE_COMMENT;
	request = std::regex_replace(request, std::regex("\\$COMMENT"), comment);
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));
	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void MySQL::delete_probe ()
{
	if (!this->probe_selected) throw std::logic_error("No point selected");

	std::string request = MySQL::DELETE_PROBE;
	request = std::regex_replace(request, std::regex("\\$POINT"), std::to_string(this->probe_id));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	this->probe_selected = false;
}

void MySQL::select_probe (double ct, double rho, double phi, double z)
{
	if (!this->problem_selected) throw std::logic_error("No problem_id selected");

	std::string request = MySQL::SELECT_PROBE;
	request = std::regex_replace(request, std::regex("\\$HEAD"), std::to_string(this->problem_id));
	request = std::regex_replace(request, std::regex("\\$TIME"), std::to_string(ct));
	request = std::regex_replace(request, std::regex("\\$RADIAL"), std::to_string(rho));
	request = std::regex_replace(request, std::regex("\\$AZIMUTH"), std::to_string(phi));
	request = std::regex_replace(request, std::regex("\\$DISTANCE"), std::to_string(z));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(this->connection);
	MYSQL_ROW row = mysql_fetch_row(responce);

	if (row != NULL) {
		this->probe_id = std::stod(row[0]);
		this->result = row[1] ? std::stod(row[1]) : NAN;
		mysql_free_result(responce);
	} else {
		this->probe_id = this->insert_probe(ct,rho,phi,z);
		this->result = NAN;
	}

	this->probe_selected = true;
}

std::size_t MySQL::insert_probe (double ct, double rho, double phi, double z)
{
	if (!this->problem_selected) throw std::logic_error("Problem is not selected");

	std::string request = MySQL::INSERT_PROBE;
	request = std::regex_replace(request, std::regex("\\$HEAD"), std::to_string(this->problem_id));
	request = std::regex_replace(request, std::regex("\\$TIME"), std::to_string(ct));
	request = std::regex_replace(request, std::regex("\\$RADIAL"), std::to_string(rho));
	request = std::regex_replace(request, std::regex("\\$AZIMUTH"), std::to_string(phi));
	request = std::regex_replace(request, std::regex("\\$DISTANCE"), std::to_string(z));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	return mysql_insert_id(this->connection);
}

std::vector<std::size_t> MySQL::select_all_probes ()
{
	if (!this->problem_selected) throw std::logic_error("Problem is not selected");

	std::string request = MySQL::SELECT_ALL_PROBES;
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(connection);
	std::size_t point_nums = mysql_num_rows(responce);

	std::vector<std::size_t> id;
	for (std::size_t item = 0; item < point_nums; item++) {
		MYSQL_ROW serch_row = mysql_fetch_row(responce);
		id.push_back(std::stoi(serch_row[0]));
	}

	mysql_free_result(responce);
	return id;
}

double MySQL::get_probe_result () const
{
	if (!this->probe_selected) throw std::logic_error("No point selected");
	return this->result;
}

void MySQL::update_probe_result (double value)
{
	if (!this->probe_selected) throw std::logic_error("No point selected");

	std::string request = MySQL::UPDATE_PROBE_RESULT;
	request = std::regex_replace(request, std::regex("\\$POINT"), std::to_string(this->probe_id));
	request = std::regex_replace(request, std::regex("\\$VALUE"), std::to_string(value));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	this->result = value;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void MySQL::delete_coeff ()
{
	if (!this->evolution_selected) throw std::logic_error("No evolution.id selected. Needed by MySQL::delete_coeff");

	std::string request = MySQL::DELETE_COEFF;
	request = std::regex_replace(request, std::regex("\\$INDEX"), std::to_string(this->evo_id));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	this->evolution_selected = false;
}

void MySQL::select_coeff (int m, double nu)
{
	if (!this->problem_selected) throw std::logic_error("No problem.id selected. Needed by MySQL::select_coeff");
	if (!this->probe_selected) throw std::logic_error("No probe.id selected. Needed by MySQL::select_coeff");

	std::string request = MySQL::SELECT_COEFF;
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));
	request = std::regex_replace(request, std::regex("\\$PROBE"), std::to_string(this->probe_id));
	request = std::regex_replace(request, std::regex("\\$M"), std::to_string(m));
	request = std::regex_replace(request, std::regex("\\$NU"), std::to_string(nu));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(this->connection);
	MYSQL_ROW row = mysql_fetch_row(responce);

	if (row != NULL) {
		this->evo_id = std::stod(row[0]);
		this->vmh = row[1] ? std::stod(row[1]) : NAN;
		mysql_free_result(responce);
	} else {
		this->evo_id = this->insert_coeff(m, nu);
		this->vmh = NAN;
	}                                                                                                                                                                                                                                                                                   

	this->evolution_selected = true;
}

std::size_t MySQL::insert_coeff (int m, double nu)
{
	if (!this->problem_selected) throw std::logic_error("No problem.id selected. Needed by MySQL::insert_coeff");
	if (!this->probe_selected) throw std::logic_error("No probe.id selected. Needed by MySQL::insert_coeff");

	std::string request = MySQL::INSERT_COEFF;
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));
	request = std::regex_replace(request, std::regex("\\$PROBE"), std::to_string(this->probe_id));
	request = std::regex_replace(request, std::regex("\\$M"), std::to_string(m));
	request = std::regex_replace(request, std::regex("\\$NU"), std::to_string(nu));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	return mysql_insert_id(this->connection);
}

std::pair<std::vector<double>,std::vector<double>> MySQL::select_all_coeffs (int m)
{
	if (!this->problem_selected) throw std::logic_error("Problem is not selected. Needed by MySQL::select_all_coeffs");
	if (!this->probe_selected) throw std::logic_error("No probe.id selected. Needed by MySQL::select_all_coeffs");

	std::string request = MySQL::SELECT_ALL_COEFF;
	request = std::regex_replace(request, std::regex("\\$PROBLEM"), std::to_string(this->problem_id));
	request = std::regex_replace(request, std::regex("\\$PROBE"), std::to_string(this->probe_id));
	if (m == 666) request = std::regex_replace(request, std::regex("\\$M_OPT"), " ");
	else request = std::regex_replace(request, std::regex("\\$M_OPT"), " AND m = " + std::to_string(m));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* responce = mysql_store_result(connection);
	std::size_t point_nums = mysql_num_rows(responce);

	std::vector<double> nu, vmh;
	nu.reserve(point_nums);
	vmh.reserve(point_nums);
	for (std::size_t item = 0; item < point_nums; item++) {
		MYSQL_ROW serch_row = mysql_fetch_row(responce);
		nu.push_back(std::stod(serch_row[0]));
		vmh.push_back(std::stod(serch_row[1]));
	}

	mysql_free_result(responce);
	return std::make_pair(nu, vmh);
}

double MySQL::get_coeff_result () const
{
	if (!this->evolution_selected) throw std::logic_error("No evolution.id selected. Needed by MySQL::get_coeff_result");
	return this->vmh;
}

void MySQL::update_coeff_result (double value)
{
	if (!this->evolution_selected) throw std::logic_error("No point selected. Needed by MySQL::update_coeff_result");

	std::string request = MySQL::UPDATE_COEFF_VAL;
	request = std::regex_replace(request, std::regex("\\$INDEX"), std::to_string(this->evo_id));
	request = std::regex_replace(request, std::regex("\\$VALUE"), std::to_string(value));

	int error_code = mysql_query(this->connection, request.c_str());
	MySQL::throw_error_code(error_code);
	this->evo_id = value;
}
