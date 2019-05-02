//
//  script_manager.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 02.05.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "logger.hpp"

#include <string>
#include <vector>
#include <utility>

class ScriptManager {

public:

    enum Colormap {gray, parula};
    ScriptManager (const std::string &filename, Logger* global_log = NULL) : fname{filename}, log{global_log} {}
    virtual ~ScriptManager () = default;

	virtual void set_title (const std::string &text) = 0;    // default: filename
	virtual void set_colormap (const Colormap &schem) = 0;   // defualt: ScriptManager::Colormap::gray
    virtual void set_ox_label (const std::string &text) = 0; // default: "OX"
	virtual void set_oy_label (const std::string &text) = 0; // default: "OY"
	virtual void set_oz_label (const std::string &text) = 0; // default: "OZ"
	virtual void set_logscale_ox (bool status = true) = 0;   // default: false
	virtual void set_logscale_oy (bool status = true) = 0;   // default: false

	virtual void plot_2d (const std::vector<double> &arg, const std::vector<double> &fnc) const = 0;
    virtual void plot_2d (const std::vector<double> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &title) const = 0;
	virtual void plot_3d (const std::vector<std::pair<double>> &args, const std::vector<double> &fnc) const = 0;
	virtual void plot_cmap (const std::vector<std::pair<double>> &args, const std::vector<double> &fnc) const = 0;

protected:
    std::string fname {};
    Logger* log {};
};
