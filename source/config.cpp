//
//  config.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 14.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "config.hpp"

bool Config::print_progress_value = true;
bool Config::call_gnuplot_value = true;
bool Config::plot_grid_value = true;
bool Config::plot_bcage_value = true;
bool Config::display_params_value = true;

Colormap Config::plot_cmap_value = Colormap::gray;
PlotDev Config::device_value = PlotDev::x11;

std::string Config::ppath_gnuplot_value = "gnuplot";
std::string Config::ppath_gnp_value = "maxwell.gmp";

std::size_t Config::h_terms_value = 100;
std::size_t Config::fbitrate_value = 256;
std::size_t Config::thread_number_value = 4;

float Config::plane_disk_radius_value = 1;
float Config::plane_disk_magnitude_value = 1;
float Config::plane_disk_epsr_value = 1;
float Config::plane_disk_mur_value = 1;

void Config::read (std::string posix_path)
{
	std::ifstream file(posix_path);
	std::string line;

	while (std::getline(file, line)) {

		std::string new_line;

		unique_copy (line.begin(), line.end(), 
			std::back_insert_iterator<std::string>(new_line),
			[] (char a, char b) { return std::isspace(a) && std::isspace(b); 
		});

		if (new_line[0] == ' ') new_line.erase(0,1);
		std::size_t lsize = new_line.size() - 1;
		if (new_line[lsize] == ' ') new_line.erase(lsize);

		if (new_line.find("#") == std::string::npos) {
			std::size_t eq_char = new_line.find("=");
			if (eq_char != std::string::npos) {

				std::string arg = new_line.substr(eq_char+1);
				if (arg[0] == ' ') arg.erase(0,1);
				std::size_t arg_size = arg.size() - 1;
				if (arg[arg_size] == ' ') arg.erase(arg_size);
				std::string option = new_line.substr(0, eq_char);

				if (option.find("PRINT_PROGRESS") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						Config::print_progress_value = true;
					else if (arg.find("FALSE") != std::string::npos)
						Config::print_progress_value = false;
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("DISPLAY_PARAMS") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						Config::display_params_value = true;
					else if (arg.find("FALSE") != std::string::npos)
						Config::display_params_value = false;
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("CALL_GNUPLOT_BIN") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						Config::call_gnuplot_value = true;
					else if (arg.find("FALSE") != std::string::npos)
						Config::call_gnuplot_value = false;
					else throw std::logic_error("Only TRUE/FALSE allowed...");
				}

				else if (option.find("PLOT_COLOR_MAP") != std::string::npos) {
					if (arg.find("GRAYSCALE") != std::string::npos)
						Config::plot_cmap_value = Colormap::gray;
					else if (arg.find("PARULA") != std::string::npos)
						Config::plot_cmap_value = Colormap::parula;
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_GRID") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						Config::plot_grid_value = true;
					else if (arg.find("FALSE") != std::string::npos)
						Config::plot_grid_value = false;
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_BOUND_CAGE") != std::string::npos) {
					if (arg.find("TRUE") != std::string::npos)
						Config::plot_bcage_value = true;
					else if (arg.find("FALSE") != std::string::npos)
						Config::plot_bcage_value = false;
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("PLOT_DEVICE") != std::string::npos) {
					if (arg.find("X11") != std::string::npos)
						Config::device_value = PlotDev::x11;
					else if (arg.find("TERM") != std::string::npos)
						Config::device_value = PlotDev::term;
					else throw std::logic_error("Not implemented...");
				}

				else if (option.find("GNUPLOT_BIN") != std::string::npos) {
					// std::experimental::filesystem::exists(arg);
					Config::ppath_gnuplot_value = arg;
				}

				else if (option.find("SCRIPT_NAME") != std::string::npos) {
					// std::experimental::filesystem::exists(arg);
					Config::ppath_gnp_value = arg;
				}

				else if (option.find("MAGNETIC_TERM_NUM") != std::string::npos) {
					Config::h_terms_value = std::stoi(arg);
				}

				else if (option.find("FLOAT_BITRATE") != std::string::npos) {
					Config::fbitrate_value = std::stoi(arg);
				}

				else if (option.find("PDISK_RADIUS") != std::string::npos) {
					Config::plane_disk_radius_value = std::stof(arg);
				}

				else if (option.find("PDISK_MAGNITUDE") != std::string::npos) {
					Config::plane_disk_magnitude_value = std::stof(arg);
				}

				else if (option.find("PDISK_EPSR") != std::string::npos) {
					Config::plane_disk_epsr_value = std::stof(arg);
				}

				else if (option.find("PDISK_MUR") != std::string::npos) {
					Config::plane_disk_mur_value = std::stof(arg);
				}

				else if (option.find("THREAD_NUMBER") != std::string::npos) {
					Config::thread_number_value = std::stof(arg);
				}
			}
		} 
	}
	file.close();
}

bool Config::print_progress ()
{
	return Config::print_progress_value;
}

bool Config::display_params ()
{
	return Config::display_params_value;
}

bool Config::call_gnuplot ()
{
	return Config::call_gnuplot_value;
}

bool Config::plot_grid ()
{
	return Config::plot_grid_value;
}

bool Config::plot_baund_cage ()
{
	return Config::plot_bcage_value;
}

Colormap Config::plot_color_map ()
{
	return Config::plot_cmap_value;
}

PlotDev Config::plot_device ()
{
	return Config::device_value;
}

std::string Config::path_gnuplot_binary ()
{
	return Config::ppath_gnuplot_value;
}

std::string Config::gnp_script_path ()
{
	return Config::ppath_gnp_value;
}

std::size_t Config::magnetic_term_num ()
{
	return Config::h_terms_value;
}

std::size_t Config::float_bitrate ()
{
	return Config::fbitrate_value;
}

std::size_t Config::thread_number ()
{
	return Config::thread_number_value;
}

// Plane Disk Problem Global Propertes

float Config::plane_disk_radius ()
{
	return Config::plane_disk_radius_value;
}

float Config::plane_disk_magnitude ()
{
	return Config::plane_disk_magnitude_value;
}

float Config::plane_disk_epsr ()
{
	return Config::plane_disk_epsr_value;
}

float Config::plane_disk_mur ()
{
	return Config::plane_disk_mur_value;
}

