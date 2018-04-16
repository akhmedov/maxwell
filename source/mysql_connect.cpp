//
//  mysql_connect.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 24.01.18.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "mysql_connect.hpp"

const std::string MySQL::SELECT_PROBLEM_ID = "SELECT id FROM maxwell.maxwell_header WHERE radiator = '$RADIATOR_TYPE' AND component = '$FIELD_COMP' AND radius = $RADIUS_VAL AND magnitude = $MAGNITUDE_VAL AND mu_r = $MUR_VAL AND eps_r = $EPSR_VAL AND kerr_r = $KERR_VAL AND duration = 0 AND signal_type = 'rect';";
const std::string MySQL::INSERT_PROBLEM    = "INSERT INTO maxwell.maxwell_header SET radiator = '$RADIATOR_TYPE', component = '$FIELD_COMP', radius = $RADIUS_VAL, magnitude = $MAGNITUDE_VAL, mu_r = $MUR_VAL, eps_r = $EPSR_VAL, kerr_r = $KERR_VAL, duration  = 0, signal_type = 'rect';";
const std::string MySQL::SELECT_POINT      = "SELECT id,lineary,square,kerr FROM maxwell.maxwell_data WHERE head_id = $HEAD AND ct = $TIME AND rho = $RADIAL AND phi = $AZIMUTH AND z = $DISTANCE;";
const std::string MySQL::INSERT_POINT      = "INSERT INTO maxwell.maxwell_data SET head_id = $HEAD, ct = $TIME, rho = $RADIAL, phi = $AZIMUTH, z = $DISTANCE;";
const std::string MySQL::UPDATE_LINEAR     = "UPDATE maxwell.maxwell_data SET lineary = $VALUE WHERE id = $POINT;";
const std::string MySQL::UPDATE_SQUARE     = "UPDATE maxwell.maxwell_data SET square  = $VALUE WHERE id = $POINT;";
const std::string MySQL::UPDATE_KERR       = "UPDATE maxwell.maxwell_data SET kerr    = $VALUE WHERE id = $POINT;";
const std::string MySQL::USE_MAXWELL       = "USE maxwell;";
const std::string MySQL::SET_WAIT_TIMEOUT  = "SET wait_timeout=9999999;";

MySQL::MySQL(Config* gl_config)
{
	this->connection = mysql_init(NULL);
	this->global_config = gl_config;

	auto connected = mysql_real_connect(this->connection,
										this->global_config->mysql_hostname().c_str(), 
										this->global_config->mysql_username().c_str(), 
										this->global_config->mysql_password().c_str(), 
										this->global_config->mysql_database().c_str(), 
										0, NULL, 0);

	if (!connected) throw std::logic_error("Connection failed! MySQL is not connected.");

	int error_code = mysql_query(this->connection, MySQL::SET_WAIT_TIMEOUT.c_str());
	MySQL::throw_error_code(error_code);

	error_code = mysql_query(this->connection, MySQL::USE_MAXWELL.c_str());
	MySQL::throw_error_code(error_code);

	std::string problem_id = MySQL::SELECT_PROBLEM_ID;
	problem_id = std::regex_replace(problem_id, std::regex("\\$RADIATOR_TYPE"), "uni_disk");
	problem_id = std::regex_replace(problem_id, std::regex("\\$FIELD_COMP"), MySQL::to_string(this->global_config->field_component()));
	problem_id = std::regex_replace(problem_id, std::regex("\\$RADIUS_VAL"), std::to_string(this->global_config->plane_disk_radius()));
	problem_id = std::regex_replace(problem_id, std::regex("\\$MAGNITUDE_VAL"), std::to_string(this->global_config->plane_disk_magnitude()));
	problem_id = std::regex_replace(problem_id, std::regex("\\$MUR_VAL"), std::to_string(this->global_config->plane_disk_mur()));
	problem_id = std::regex_replace(problem_id, std::regex("\\$EPSR_VAL"), std::to_string(this->global_config->plane_disk_epsr()));
	problem_id = std::regex_replace(problem_id, std::regex("\\$KERR_VAL"), std::to_string(this->global_config->kerr_value()));

	error_code = mysql_query(this->connection, problem_id.c_str());
	MySQL::throw_error_code(error_code);
	MYSQL_RES* result_init = mysql_store_result(connection);
	MYSQL_ROW row_init = mysql_fetch_row(result_init);

	if (row_init != NULL) {

		this->problem_id = std::stoi(row_init[0]);
		mysql_free_result(result_init);

	} else {

		std::string insert_problem = MySQL::INSERT_PROBLEM;
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$RADIATOR_TYPE"), "uni_disk");
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$FIELD_COMP"), MySQL::to_string(this->global_config->field_component()));
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$RADIUS_VAL"), std::to_string(this->global_config->plane_disk_radius()));
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$MAGNITUDE_VAL"), std::to_string(this->global_config->plane_disk_magnitude()));
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$MUR_VAL"), std::to_string(this->global_config->plane_disk_mur()));
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$EPSR_VAL"), std::to_string(this->global_config->plane_disk_epsr()));
		insert_problem = std::regex_replace(insert_problem, std::regex("\\$KERR_VAL"), std::to_string(this->global_config->kerr_value()));

		error_code = mysql_query(this->connection, insert_problem.c_str());
		MySQL::throw_error_code(error_code);

		error_code = mysql_query(this->connection, problem_id.c_str());
		MySQL::throw_error_code(error_code);

		mysql_free_result(result_init);
		MYSQL_RES* result_new = mysql_store_result(this->connection);
		MYSQL_ROW row_new = mysql_fetch_row(result_new);

		if (row_new != NULL) {
			this->problem_id = std::stoi(row_new[0]);
			mysql_free_result(result_new);
		} else throw std::logic_error("Internal maxwell error in MySQL module");
	}
}

void MySQL::select_point (double ct, double rho, double phi, double z)
{
	std::string select_point = MySQL::SELECT_POINT;
	select_point = std::regex_replace(select_point, std::regex("\\$HEAD"), std::to_string(this->problem_id));
	select_point = std::regex_replace(select_point, std::regex("\\$TIME"), std::to_string(ct));
	select_point = std::regex_replace(select_point, std::regex("\\$RADIAL"), std::to_string(rho));
	select_point = std::regex_replace(select_point, std::regex("\\$AZIMUTH"), std::to_string(phi));
	select_point = std::regex_replace(select_point, std::regex("\\$DISTANCE"), std::to_string(z));

	int error_code = mysql_query(this->connection, select_point.c_str());
	MySQL::throw_error_code(error_code);

	MYSQL_RES* serch_result = mysql_store_result(this->connection);
	MYSQL_ROW serch_row = mysql_fetch_row(serch_result);

	if (serch_row != NULL) {

		this->point_id = std::stod(serch_row[0]);
		this->linear   = serch_row[1] ? std::stod(serch_row[1]) : NAN;
		this->square   = serch_row[2] ? std::stod(serch_row[2]) : NAN;
		this->kerr     = serch_row[3] ? std::stod(serch_row[3]) : NAN;
		mysql_free_result(serch_result);

	} else {

		std::string insert_point = MySQL::INSERT_POINT;
		insert_point = std::regex_replace(insert_point, std::regex("\\$HEAD"), std::to_string(this->problem_id));
		insert_point = std::regex_replace(insert_point, std::regex("\\$TIME"), std::to_string(ct));
		insert_point = std::regex_replace(insert_point, std::regex("\\$RADIAL"), std::to_string(rho));
		insert_point = std::regex_replace(insert_point, std::regex("\\$AZIMUTH"), std::to_string(phi));
		insert_point = std::regex_replace(insert_point, std::regex("\\$DISTANCE"), std::to_string(z));

		error_code = mysql_query(this->connection, insert_point.c_str());
		MySQL::throw_error_code(error_code);
		mysql_free_result(serch_result);

		/* error_code = mysql_query(this->connection, select_point.c_str());
		MySQL::throw_error_code(error_code);
		serch_result = mysql_store_result(this->connection);
		serch_row = mysql_fetch_row(serch_result);

		if (serch_row != NULL) {
			this->point_id = std::stod(serch_row[0]);
			this->linear   = serch_row[1] ? std::stod(serch_row[1]) : NAN;
			this->square   = serch_row[2] ? std::stod(serch_row[2]) : NAN;
			this->kerr     = serch_row[3] ? std::stod(serch_row[3]) : NAN;
			mysql_free_result(serch_result);
		} else throw std::logic_error("Internal maxwell error in MySQL module"); */

		/* LAST_INSERTED_ID() with multi threading and several connected cliets?
		The ID that was generated is maintained in the server on a per-connection 
		basis. This means that the value returned by the function to a given 
		client is the first AUTO_INCREMENT value generated for most recent 
		statement affecting an AUTO_INCREMENT column by that client. This value 
		cannot be affected by other clients, even if they generate AUTO_INCREMENT 
		values of their own. This behavior ensures that each client can retrieve 
		its own ID without concern for the activity of other clients, and without 
		the need for locks or transactions. So unless your inserts for multiple 
		users would happen to be made over the same database connection, you have 
		nothing to worry about. */

		this->point_id = mysql_insert_id(this->connection);
		this->linear   = NAN;
		this->square   = NAN;
		this->kerr     = NAN;
	}
}

std::string MySQL::get_hostname() const
{
	return this->global_config->mysql_username() + "@" + this->global_config->mysql_hostname();
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
		default: message += std::to_string(code);
	}
	throw std::logic_error(message);
}

std::string MySQL::to_string(const FieldComponent& type)
{
	switch (type) {
		case 0: return "Ex";
		case 1: return "Ey";
		case 2: return "Ez";
		case 3: return "Ephi";
		case 4: return "Erho";
		case 5: return "Hx";
		case 6: return "Hy";
		case 7: return "Hz";
		case 8: return "Hphi";
		case 9: return "Hrho";
	}
}

MySQL::~MySQL ()
{
	mysql_close(this->connection);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	delete this->global_config;
	delete this->connection;
}

//=======================================================================================

double MySQL::get_linear () const
{
	return this->linear;
}

void MySQL::set_linear (double value)
{
	this->linear = value;
	
	std::string update_linear = MySQL::UPDATE_LINEAR;
	update_linear = std::regex_replace(update_linear, std::regex("\\$VALUE"), std::to_string(value));
	update_linear = std::regex_replace(update_linear, std::regex("\\$POINT"), std::to_string(this->point_id));
	
	int error_code = mysql_query(this->connection, update_linear.c_str());
	MySQL::throw_error_code(error_code);
}

//=======================================================================================

double MySQL::get_square () const
{
	return this->square;
}

void MySQL::set_square (double value)
{
	this->square = value;

	std::string update_square = MySQL::UPDATE_SQUARE;
	update_square = std::regex_replace(update_square, std::regex("\\$VALUE"), std::to_string(value));
	update_square = std::regex_replace(update_square, std::regex("\\$POINT"), std::to_string(this->point_id));

	int error_code = mysql_query(this->connection, update_square.c_str());
	MySQL::throw_error_code(error_code);
}

//=======================================================================================

double MySQL::get_kerr () const
{
	return this->kerr;
}

void MySQL::set_kerr (double value)
{
	this->kerr = value;
	
	std::string update_kerr = MySQL::UPDATE_KERR;
	update_kerr = std::regex_replace(update_kerr, std::regex("\\$VALUE"), std::to_string(value));
	update_kerr = std::regex_replace(update_kerr, std::regex("\\$POINT"), std::to_string(this->point_id));
	
	int error_code = mysql_query(this->connection, update_kerr.c_str());
	MySQL::throw_error_code(error_code);
}

//=======================================================================================

double MySQL::get_value(const std::type_info& type) const
{
	if (type == typeid(MissileField)) {
		return this->get_linear();
	} else if (type == typeid(KerrAmendment)) {
		return this->get_kerr();
	} else {
		std::string type_name = type.name();
		throw std::logic_error(type_name + " is not implemented in MySQL::get_value");
	}
}

void MySQL::set_value(const std::type_info& type, double value)
{
	if (type == typeid(MissileField)) {
		this->set_linear(value);
		return;
	} else if (type == typeid(KerrAmendment)) {
		this->set_kerr(value);
		return;
	} else {
		std::string type_name = type.name();
		throw std::logic_error(type_name + " is not implemented in MySQL::set_value");
	}
}
