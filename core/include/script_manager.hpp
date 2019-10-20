//
//  script_manager.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 02.05.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "logger.hpp"

#include <map>
#include <regex>
#include <string>
#include <vector>
#include <utility>
#include <fstream>

class ScriptManager {

public:

	enum Colormap {grey, jet, hot, coolwarm};
	ScriptManager (const std::string &filename, Logger* global_log = NULL) : fname{filename}, log{global_log} {}
	virtual ~ScriptManager () = default;

	void set_colormap (const Colormap &schem) { this->color = schem; }
	void set_title (const std::string &text) { this->title = text; }
	void set_fontsize (size_t val) { this->font = std::to_string(val); }
	void set_ox_label (const std::string &text) { this->ox = text; }
	void set_oy_label (const std::string &text) { this->oy = text; }
	void set_oz_label (const std::string &text) { this->oz = text; }
	void set_logscale_ox (bool status) { this->logX = status; }
	void set_logscale_oy (bool status) { this->logY = status; }

	virtual void plot2d (const std::vector<double> &arg, const std::vector<double> &fnc) const = 0;
	virtual void plot2d (const std::vector<std::vector<double>> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &title) const = 0;
	virtual void plot3d (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const = 0;
	virtual void colormap (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const = 0;

protected:

	static std::vector<std::vector<double>> datagrid_from (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc, double &xmin, double &xstep, double &xmax, double &ymin, double &ystep, double &ymax);
	void write_script (const std::string &text) const;

	std::string fname {"script.py"};
	Logger* log {NULL};

	Colormap color {Colormap::grey};
	std::string font {"14"};
	std::string title {""};
	std::string ox {"OX"};
	std::string oy {"OY"};
	std::string oz {"OZ"};
	bool logX {false};
	bool logY {false};
};
