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
	text_file.open(this->file_name);
	if (!this->text_file) throw std::logic_error("Logging file is not open.");
	this->info("Logging sistem is started.");
}

void Logger::info (const std::string &text)
{
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::II, text);
	query.push(mesg);
}

void Logger::warning (const std::string &text)
{
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::WW, text);
	query.push(mesg);
}

void Logger::error (const std::string &text)
{
	if (this->query.size() >= QSIZE_LIMIT) this->write_all();
	auto mesg = std::make_pair(Logger::MesgType::EE, text);
	query.push(mesg);
}

void Logger::write_next ()
{
	if (!this->text_file) {
		this->error("Log file unavailable. Reopenning...");
		text_file.open(this->file_name);
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
