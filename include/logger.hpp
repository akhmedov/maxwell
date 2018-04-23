//
//  logger.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 20.04.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef logger_hpp
#define logger_hpp

#define QSIZE_LIMIT 1

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#include <mutex>
#include <ctime>
#include <queue>
#include <string>
#include <utility>
#include <fstream>

struct Logger {
	Logger (const std::string &posix_path);
	void info (const std::string &text);
	void warning (const std::string &text);
	void error (const std::string &text);
	void write_next ();
	void write_all ();
	~Logger ();

protected:
	static std::string sys_time ();
	static std::string name_upgrade (const std::string &old_name);

private:
	std::mutex mesg_add;
	const std::string file_name;
	enum MesgType {II, WW, EE};
	std::ofstream text_file;
	std::queue<std::pair<MesgType,std::string>> query;
};

#endif /* logger_hpp */
