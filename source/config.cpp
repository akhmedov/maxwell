//
//  config.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 14.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "config.hpp"

Config::Config ()
{
	this->plot_cl_option = 0;
	this->dataset_cl_option = 0;
	this->version_cl_option = false;
	this->help_cl_option = false;

	this->print_progress_value = true;
	this->call_gnuplot_value = true;
	this->plot_grid_value = true;
	this->plot_bcage_value = true;
	this->display_params_value = true;
	this->is_safe_mode = false;
	this->is_logger_mode = false;

	this->plot_cmap_value = Colormap::gray;
	this->device_value = PlotDev::x11;

	this->ppath_gnuplot_value = "gnuplot";
	this->ppath_gnp_value = "maxwell.gnp";
	this->ppath_maxwell_config = "maxwell.conf";
	this->ppath_maxwell_log = "maxwell.log";

	this->h_terms_value = 100;
	this->fbitrate_value = 256;
	this->thread_number_value = 4;

	this->impulse_shape_value = ImpulseShape::rect;

	this->plane_disk_radius_value = 1;
	this->plane_disk_magnitude_value = 1;
	this->plane_disk_epsr_value = 1;
	this->plane_disk_mur_value = 1;
	this->kerr_coeff_value = 0;
	this->is_medium_kerr = false;
	this->noise_percent = 0;
	this->superposition = Superposition::additive;
	
	// working_component must not be default!!!
	// this->working_component = FieldComponnt::Ex;

	this->mysql_addr  = "localhost";
	this->mysql_user  = "maxwell";
	this->mysql_pass  = "maxwell";
	this->mysql_dbase = "maxwell";

	/* TODO: Non default options
	this->vt_value;
	this->rho_value;
	this->phi_value;
	this->z_value; */
}

/* setters */

void Config::receiver_vt (double value)
{
	this->vt_value = {{value, value, value}};
}

void Config::receiver_rho (double value)
{
	this->rho_value = {{value, value, value}};
}

void Config::receiver_phi (double value)
{
	this->phi_value = {{value, value, value}};
}

void Config::receiver_z (double value)
{
	this->z_value = {{value, value, value}};
}

void Config::receiver_vt (double from, double step, double to)
{
	this->vt_value = {{from, step, to}};
}

void Config::receiver_rho (double from, double step, double to)
{
	this->rho_value = {{from, step, to}};
}

void Config::receiver_phi (double from, double step, double to)
{
	this->phi_value = {{from, step, to}};
}

void Config::receiver_z (double from, double step, double to)
{
	this->z_value = {{from, step, to}};
}

void Config::version_opt (bool option)
{
	this->version_cl_option = option;
}

void Config::help_opt (bool option)
{
	this->help_cl_option = option;
}

void Config::dataset_model (std::size_t enum_model_name)
{
	this->dataset_cl_option = enum_model_name;
}

void Config::plot_model (std::size_t enum_model_name)
{
	this->plot_cl_option = enum_model_name;
}

void Config::print_progress (bool option)
{
	this->print_progress_value = option;
}

void Config::plot_grid (bool option)
{
	this->plot_grid_value = option;
}

void Config::plot_baund_cage (bool option)
{
	this->plot_bcage_value = option;
}

void Config::call_gnuplot (bool option) 
{
	this->call_gnuplot_value = option;
}

void Config::display_params (bool option)
{
	this->display_params_value = option;
}

void Config::safe_mode (bool option)
{
	this->is_safe_mode = option;
}

void Config::logger_status (bool option)
{
	this->is_logger_mode = option;
}

void Config::plot_color_map (Colormap plot_color)
{
	this->plot_cmap_value = plot_color;
}

void Config::plot_device (PlotDev plot_device)
{
	this->device_value = plot_device;
}

void Config::path_gnuplot_binary (std::string posix_path)
{
	this->ppath_gnuplot_value = posix_path;
}

void Config::gnp_script_path (std::string posix_path)
{
	this->ppath_gnp_value = posix_path;
}

void Config::maxwell_config_path (std::string posix_path)
{
	this->ppath_maxwell_config = posix_path;
}

void Config::maxwell_log_path (std::string posix_path)
{
	this->ppath_maxwell_log = posix_path;
}

void Config::magnetic_term_num (std::size_t terms)
{
	this->h_terms_value = terms;
}

void Config::float_bitrate (std::size_t float_digits)
{
	this->fbitrate_value = float_digits;
}

void Config::thread_number (std::size_t threads)
{
	this->thread_number_value = threads;
}

void Config::impulse_shape (ImpulseShape shape)
{
	this->impulse_shape_value = shape;
}

void Config::plane_disk_radius (float radius)
{
	this->plane_disk_radius_value = radius;
}

void Config::plane_disk_magnitude (float magnitude)
{
	this->plane_disk_magnitude_value = magnitude;
}

void Config::plane_disk_epsr (float eps_r)
{
	this->plane_disk_epsr_value = eps_r;
}

void Config::plane_disk_mur (float mu_r)
{
	this->plane_disk_mur_value = mu_r;
}

void Config::kerr_value (double chi_3)
{
	this->kerr_coeff_value = chi_3;
}

void Config::kerr_medium (bool option)
{
	this->is_medium_kerr = option;
}

void Config::noise_level (double persent)
{
	this->noise_percent = persent;
}

void Config::medium_superposition (Superposition type)
{
	this->superposition = type;
}

void Config::field_component (std::size_t model_num)
{
	switch (model_num) {
		case 1: this->working_component = FieldComponent::Ex; break;
		case 2: this->working_component = FieldComponent::Hy; break;
		case 3: this->working_component = FieldComponent::Ex; break;
		case 4: this->working_component = FieldComponent::Hy; break;
		case 5: this->working_component = FieldComponent::Ex; break;
		case 6: this->working_component = FieldComponent::Hy; break;
		case 7: this->working_component = FieldComponent::Hy; break;
		default: throw std::logic_error("This model number is not implemnted in Config::field_component()");

	};
}

void Config::mysql_hostname (std::string text)
{
	this->mysql_addr = text;
}

void Config::mysql_username (std::string text)
{
	this->mysql_user = text;
}

void Config::mysql_password (std::string text)
{
	this->mysql_pass = text;
}

void Config::mysql_database (std::string text)
{
	this->mysql_dbase = text;
}

/* getters */

double Config::kerr_value () const
{
	return this->kerr_coeff_value;
}

bool Config::kerr_medium () const
{
	return this->is_medium_kerr;
}

std::array<double,3> Config::receiver_vt () const
{
	return this->vt_value;
}

std::array<double,3> Config::receiver_rho () const
{
	return this->rho_value;
}

std::array<double,3> Config::receiver_phi () const
{
	return this->phi_value;
}

std::array<double,3> Config::receiver_z () const
{
	return this->z_value;
}

bool Config::version_opt () const
{
	return this->version_cl_option; 
}

bool Config::help_opt () const
{
	return this->help_cl_option;
}

std::size_t Config::dataset_model () const
{
	return this->dataset_cl_option;
}

std::size_t Config::plot_model () const
{
	return this->plot_cl_option;
}

bool Config::print_progress () const
{
	return this->print_progress_value;
}

bool Config::display_params () const
{
	return this->display_params_value;
}

bool Config::call_gnuplot () const
{
	return this->call_gnuplot_value;
}

bool Config::plot_grid () const
{
	return this->plot_grid_value;
}

bool Config::plot_baund_cage () const
{
	return this->plot_bcage_value;
}

bool Config::safe_mode () const
{
	return this->is_safe_mode;
}

bool Config::logger_status () const
{
	return this->is_logger_mode;
}

Colormap Config::plot_color_map () const
{
	return this->plot_cmap_value;
}

PlotDev Config::plot_device () const
{
	return this->device_value;
}

std::string Config::path_gnuplot_binary () const
{
	return this->ppath_gnuplot_value;
}

std::string Config::gnp_script_path () const
{
	return this->ppath_gnp_value;
}

std::string Config::maxwell_config_path () const
{
	return this->ppath_maxwell_config;
}

std::string Config::maxwell_log_path () const
{
	return this->ppath_maxwell_log;
}

std::size_t Config::magnetic_term_num () const
{
	return this->h_terms_value;
}

std::size_t Config::float_bitrate () const
{
	return this->fbitrate_value;
}

std::size_t Config::thread_number () const
{
	return this->thread_number_value;
}

// Plane Disk Problem Global Propertes

ImpulseShape Config::impulse_shape () const
{
	return this->impulse_shape_value;
}

float Config::plane_disk_radius () const
{
	return this->plane_disk_radius_value;
}

float Config::plane_disk_magnitude () const
{
	return this->plane_disk_magnitude_value;
}

float Config::plane_disk_epsr () const
{
	return this->plane_disk_epsr_value;
}

float Config::plane_disk_mur () const
{
	return this->plane_disk_mur_value;
}

double Config::noise_level () const
{
	return this->noise_percent;
}

Superposition Config::medium_superposition () const
{
	return this->superposition;
}

FieldComponent Config::field_component () const
{
	return this->working_component;
}

std::string Config::mysql_hostname () const
{
	return this->mysql_addr;
}

std::string Config::mysql_username () const
{
	return this->mysql_user;
}

std::string Config::mysql_password () const
{
	return this->mysql_pass;
}

std::string Config::mysql_database () const
{
	return this->mysql_dbase;
}
