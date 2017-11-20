//
//  config.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef config_hpp
#define config_hpp

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

enum Colormap {gray, parula};
enum PlotDev {x11, term};

struct Config {

	static void read (std::string posix_path);

	static bool print_progress ();
	static bool plot_grid ();
	static bool plot_baund_cage ();
	static bool call_gnuplot ();
	static bool display_params ();

	static Colormap plot_color_map ();
	static PlotDev plot_device ();

	static std::string path_gnuplot_binary ();
	static std::string gnp_script_path ();

	static std::size_t magnetic_term_num ();
	static std::size_t float_bitrate ();
	static std::size_t thread_number ();

	static float plane_disk_radius ();
	static float plane_disk_magnitude ();
	static float plane_disk_epsr ();
	static float plane_disk_mur ();

private:
	static bool print_progress_value;
	static bool plot_grid_value;
	static bool plot_bcage_value;
	static bool call_gnuplot_value;
	static bool display_params_value;

	static Colormap plot_cmap_value;
	static PlotDev device_value;

	static std::string ppath_gnuplot_value;
	static std::string ppath_gnp_value;

	static std::size_t h_terms_value;
	static std::size_t fbitrate_value;
	static std::size_t thread_number_value;

	static float plane_disk_radius_value;
	static float plane_disk_magnitude_value;
	static float plane_disk_epsr_value;
	static float plane_disk_mur_value;
};

#endif /* config_hpp */
