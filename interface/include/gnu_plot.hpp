//
//  gnu_plot.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 30.06.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef gnu_plot_hpp
#define gnu_plot_hpp

#include <cmath>
#include <regex>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <cstdlib>
#include <fstream>
#include <cstddef>
#include <string>
#include <vector>
#include <list>
#include <map>

#include "config.hpp"

typedef std::vector<std::string> Text;

template <typename T> 
using Plot2D = std::vector<std::pair<T, T>>;

template <typename T> 
using PlotMulti = std::vector<std::vector<T>>;

template <typename T> 
using Plot3D = std::vector<std::tuple<T, T, T>>;

/*
-- Data structure for GnuPlot:: module 
-- plot2d     		data: x f
-- plot3d     		data: x y f
-- plot_multi 		data: x y f1 f2 f3 ...
-- plot_colormap	data: {{f(x1,y1), f(x1,y2)}, {f(x2,y1), f(x2,y2)}}
*/

struct GnuPlot {
	GnuPlot (std::string filename);
	void set_gnuplot_bin (std::string filename);
	void set_font_size (std::size_t size);
	void set_ox_label (std::string text);
	void set_oy_label (std::string text);
	void set_oz_label (std::string text);
	void set_title (std::string text);
	void set_colormap (const Colormap &schem);
	void set_logscale_ox (bool status = true);
	void set_logscale_oy (bool status = true);
	void grid_on (bool status = true);
	void cage_on (bool status = true);

	void plot2d (const std::vector<std::vector<double>> &array);
	void plot_multi (const std::vector<std::vector<double>> &array, 
					 const std::vector<std::string> &title);
	void plot3d (const std::vector<std::vector<double>> &matrix);
	void plot_colormap (const std::vector<std::vector<double>> &array, int axis1, int axis2);
	void call_gnuplot ();

protected:

	static std::vector<std::vector<double>> datagrid_from (const std::vector<std::vector<double>>& array, int axis1, int axis2);

private:
	void direct_gnuplot_call (const Text &plot_data) const;
	void write_script (const std::string &plot_data) const;
	void write_commants_to_script (const Text &plot_data, const std::vector<std::string> &title = {}) const;

	// templates

	const static std::string MULTILINE;
	const static std::string LINE;
	const static std::string SURFACE;
	const static std::string CMAP;

	// variables 

	std::string gnuplot_path;
	std::string script_name;

	std::string font;
	std::string title;
	std::string ox_label;
	std::string oy_label;
	std::string oz_label;

	Colormap color_schem;

	bool is_multi_plot;
	bool is_title_on;
	bool is_ox_logscale, is_oy_logscale;
	bool is_grid_on;
	bool is_plot_in_cage;
	bool is_3d_plot;
	bool is_colormap_plot;
};

#endif /* gnu_plot_hpp */
