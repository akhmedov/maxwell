//
//  gnuplot_manager.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 30.06.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef gnuplot_manager_hpp
#define gnuplot_manager_hpp

#include "script_manager.hpp"

#include <regex>
#include <string>
#include <fstream>
#include <vector>
#include <map>

struct GnuPlotManager : public ScriptManager {

	GnuPlotManager (const std::string &filename, Logger* global_log = NULL);
	void plot2d (const std::vector<double> &arg, const std::vector<double> &fnc) const;
	void plot2d (const std::vector<std::vector<double>> &arg, const std::vector<std::vector<double>> &fnc, const std::vector<std::string> &title) const;
	void plot3d (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const;
	void colormap (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc) const;

private:

	const static std::string MULTILINE;
	const static std::string LINE;
	const static std::string SURFACE;
	const static std::string CMAP;
};

#endif /* gnuplot_manager_hpp */
