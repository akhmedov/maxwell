//
//  pyplot_manager.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 15.10.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "script_manager.hpp"

#include <cmath>
#include <string>
#include <vector>
#include <utility>

class PyPlotManager : public ScriptManager {

public:

	PyPlotManager (const std::string &filename, Logger* global_log = NULL);
	void plot2d (const std::vector<double> &arg, const std::vector<double> &fnc) const;
	void plot2d (const std::vector<std::vector<double>> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &title) const;
	void plot3d (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const;
	void colormap (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const;

private:

	const static std::string MULTILINE;
	const static std::string SURFACE;
	const static std::string CMAP;

};
