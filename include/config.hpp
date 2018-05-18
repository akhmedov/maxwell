//
//  config.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 13.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef config_hpp
#define config_hpp

#include <string>
#include <array>

enum Superposition {additive, multipl};
enum ModelType {plot, dataset, info};
enum ImpulseShape {on, meandr, duhamel};
enum Colormap {gray, parula};
enum PlotDev {x11, term};
enum FieldComponent {Ex, Ey, Ez, Ephi, Erho, Hx, Hy, Hz, Hphi, Hrho};

struct Config {

	Config ();

	/* setters */

	void receiver_vt (double value);
	void receiver_rho (double value);
	void receiver_phi (double value);
	void receiver_z (double value);

	void receiver_vt (double from, double step, double to);
	void receiver_rho (double from, double step, double to);
	void receiver_phi (double from, double step, double to);
	void receiver_z (double from, double step, double to);

	void version_opt (bool);
	void help_opt (bool);
	void dataset_model (std::size_t);
	void plot_model (std::size_t);

	void print_progress (bool);
	void plot_grid (bool);
	void plot_baund_cage (bool);
	void call_gnuplot (bool);
	void display_params (bool);
	void reset_noise (bool);
	void safe_mode (bool);
	void logger_status (bool);

	void plot_color_map (Colormap);
	void plot_device (PlotDev);

	void path_gnuplot_binary (std::string);
	void gnp_script_path (std::string);
	void maxwell_config_path (std::string);
	void maxwell_log_path (std::string);

	void magnetic_term_num (std::size_t);
	void float_bitrate (std::size_t);
	void thread_number (std::size_t);

	void impulse_shape (ImpulseShape);
	void plane_disk_radius (float);
	void plane_disk_magnitude (float);
	void plane_disk_epsr (float);
	void plane_disk_mur (float);
	void kerr_value (double);
	void kerr_medium (bool);
	void noise_level (double);
	void medium_superposition (Superposition);
	void field_component (std::size_t model_num);
	void duration (double);

	void mysql_hostname (std::string);
	void mysql_username (std::string);
	void mysql_password (std::string);
	void mysql_database (std::string);

	/* getters */

	std::array<double,3> receiver_vt () const;
	std::array<double,3> receiver_rho () const;
	std::array<double,3> receiver_phi () const;
	std::array<double,3> receiver_z () const;

	bool version_opt () const;
	bool help_opt () const;
	std::size_t dataset_model () const;
	std::size_t plot_model () const;

	bool print_progress () const;
	bool plot_grid () const;
	bool plot_baund_cage () const;
	bool call_gnuplot () const;
	bool display_params () const;
	bool reset_noise () const;
	bool safe_mode () const;
	bool logger_status () const;

	Colormap plot_color_map () const;
	PlotDev plot_device () const;

	std::string path_gnuplot_binary () const;
	std::string gnp_script_path () const;
	std::string maxwell_config_path () const;
	std::string maxwell_log_path () const;

	std::size_t magnetic_term_num () const;
	std::size_t float_bitrate () const;
	std::size_t thread_number () const;

	ImpulseShape impulse_shape () const;
	float plane_disk_radius () const;
	float plane_disk_magnitude () const;
	float plane_disk_epsr () const;
	float plane_disk_mur () const;
	double kerr_value () const;
	bool kerr_medium () const;
	double noise_level () const;
	Superposition medium_superposition () const;
	FieldComponent field_component () const;
	double duration () const;

	std::string mysql_hostname () const;
	std::string mysql_username () const;
	std::string mysql_password () const;
	std::string mysql_database () const;

private:

	// TODO: update with c++17 std::variant
	// TODO: receiver world point values
	std::array<double,3> vt_value;
	std::array<double,3> rho_value;
	std::array<double,3> phi_value;
	std::array<double,3> z_value;

	bool version_cl_option;
	bool help_cl_option;
	std::size_t dataset_cl_option;
	std::size_t plot_cl_option;

	bool print_progress_value;
	bool plot_grid_value;
	bool plot_bcage_value;
	bool call_gnuplot_value;
	bool display_params_value;
	bool is_safe_mode;
	bool is_logger_mode;

	Colormap plot_cmap_value;
	PlotDev device_value;

	std::string ppath_gnuplot_value;
	std::string ppath_gnp_value;
	std::string ppath_maxwell_config;
	std::string ppath_maxwell_log;

	std::size_t h_terms_value;
	std::size_t fbitrate_value;
	std::size_t thread_number_value;

	ImpulseShape impulse_shape_value;
	float plane_disk_radius_value;
	float plane_disk_magnitude_value;
	float plane_disk_epsr_value;
	float plane_disk_mur_value;
	double kerr_coeff_value;
	bool is_medium_kerr;
	double noise_percent;
	Superposition superposition;
	FieldComponent working_component;
	double signal_duration;

	std::string mysql_addr;
	std::string mysql_user;
	std::string mysql_pass;
	std::string mysql_dbase;
};

#endif /* config_hpp */
