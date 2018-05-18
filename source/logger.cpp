//
//  logger.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 20.04.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include"logger.hpp"

Logger::Logger (const std::string &posix_path)
: file_name(posix_path) {
	text_file.open(Logger::name_upgrade(posix_path));
	if (!this->text_file) throw std::logic_error("Logging file is not open.");
	this->info("Logging sistem is started.");
}

void Logger::info (const std::string &text)
{
	this->mesg_add.lock();
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::II, text);
	query.push(mesg);
	this->mesg_add.unlock();
}

void Logger::warning (const std::string &text)
{
	this->mesg_add.lock();
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::WW, text);
	query.push(mesg);
	this->mesg_add.unlock();
}

void Logger::error (const std::string &text)
{
	this->mesg_add.lock();
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::EE, text);
	query.push(mesg);
	this->mesg_add.unlock();
}

void Logger::write_next ()
{
	if (!this->text_file) {
		this->text_file.close();
		this->error("Log file unavailable. Reopenning...");
		text_file.open(Logger::name_upgrade(this->file_name));
	}

	auto mesg = this->query.front();
	this->query.pop();

	switch (mesg.first) {
		case Logger::MesgType::II: this->text_file << "[II]"; break;
		case Logger::MesgType::WW: this->text_file << "[WW]"; break;
		case Logger::MesgType::EE: this->text_file << "[EE]"; break;
		default: throw std::logic_error("Undefined Logger::MesgType catched.");
	}

	this->text_file << ' ' << mesg.second << std::endl;
}

void Logger::write_all ()
{
	while (!this->query.empty()) this->write_next();
}

std::string Logger::sys_time ()
{
	std::time_t rawtime;
	std::tm* timeinfo;
	char buffer [80];

	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);

	std::strftime(buffer,80,"%Y%m%d-%H%M%S",timeinfo);
	// std::puts(buffer);
	return std::string(buffer);
}

std::string Logger::name_upgrade (const std::string &basic)
{
	const std::string extention = ".log";
	std::size_t end_pos = basic.find(extention);

	if (end_pos != std::string::npos) {
		std::string res = basic;
		std::string add = "-" + Logger::sys_time();
		res.insert(end_pos, add);
		return res;
	}

	throw std::logic_error("Illegal logging file format. Must ends with .log extention.");
}

Logger::~Logger ()
{
	this->info("Logger destructor is called");
	this->write_all();
	text_file.close();
}
