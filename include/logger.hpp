//
//  logger.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 20.04.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef logger_hpp
#define logger_hpp

#define QSIZE_LIMIT 5e3

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#include <string>
#include <queue>
#include <utility>
#include <fstream>

struct Logger {
	Logger (const std::string &posix_path);
	void info (const std::string &text);
	void warning (const std::string &text);
	void error (const std::string &text);
	void write_next ();
	void write_all ();
private:
	const std::string file_name;
	enum MesgType {II, WW, EE};
	std::ofstream text_file;
	std::queue<std::pair<MesgType,std::string>> query;
};

#endif /* logger_hpp */
