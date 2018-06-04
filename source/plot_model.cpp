//
//  plot_model.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 27.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "plot_model.hpp"

PlotModel::PlotModel (Config* conf)
{
	this->global_log = NULL;
	this->global_conf = conf;
	this->model_pointer = {
		&PlotModel::__from_ct,
		&PlotModel::__from_ct_rho,
		&PlotModel::__from_ct_z,
		&PlotModel::__from_x_y
	};
}

void PlotModel::set_logger (Logger* log)
{
	this->global_log = log;
}

void PlotModel::call (const PlotModel::Name& model_name, const std::vector<std::pair<Component,AbstractField*>>& to_compute)
{
	this->model_pointer[model_name - 1](*this, to_compute);
}

void PlotModel::__from_ct (const std::vector<std::pair<Component,AbstractField*>>& to_compute)
{
	std::size_t thread_num = this->global_conf->thread_number();
	double ct_from = this->global_conf->receiver_vt()[0];
	double ct_step = this->global_conf->receiver_vt()[1];
	double ct_to = this->global_conf->receiver_vt()[2];
	double rho = this->global_conf->receiver_rho()[0];
	double phi = this->global_conf->receiver_phi()[0];
	double z = this->global_conf->receiver_z()[0];

	Manager* thead_core;
	if (!this->global_conf->safe_mode()) thead_core = new Manager(thread_num, this->global_log);
	else thead_core = new SafeManager( thread_num, this->global_conf, this->global_log);
	thead_core->progress_bar( this->global_conf->print_progress() );

	for (double ct = ct_from; ct <= ct_to; ct += ct_step)
		thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call(to_compute);
	std::vector<std::vector<double>> data = thead_core->get_value();

	// Superposition medium_type = this->global_conf->medium_superposition();
	std::vector<std::vector<double>> plot_data;
	for (auto i : data) {
		double magnitude = 0;
		for (std::size_t j = 4; j < i.size(); j++) magnitude += i[j];
		plot_data.push_back({i[0], magnitude});
	}

	GnuPlot* plot = new GnuPlot(this->global_conf->gnp_script_path());
	plot->set_gnuplot_bin(this->global_conf->path_gnuplot_binary());
	plot->set_colormap(this->global_conf->plot_color_map());
	plot->set_ox_label("ct, m");
	plot->set_oy_label("Ex, V/m");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(plot_data);
	if ( this->global_conf->call_gnuplot() ) plot->call_gnuplot();
}

void PlotModel::__from_ct_rho (const std::vector<std::pair<Component,AbstractField*>>& to_compute)
{
	std::size_t thread_num = this->global_conf->thread_number();
	double ct_from = this->global_conf->receiver_vt()[0];
	double ct_step = this->global_conf->receiver_vt()[1];
	double ct_to = this->global_conf->receiver_vt()[2];
	double rho_from = this->global_conf->receiver_rho()[0];
	double rho_step = this->global_conf->receiver_rho()[1];
	double rho_to = this->global_conf->receiver_rho()[2];
	double phi = this->global_conf->receiver_phi()[0];
	double z = this->global_conf->receiver_z()[0];

	Manager* thead_core;
	if (!this->global_conf->safe_mode()) thead_core = new Manager( thread_num );
	else thead_core = new SafeManager( thread_num, this->global_conf );
	thead_core->progress_bar( this->global_conf->print_progress() );

	for (double rho = rho_from; rho <= rho_to; rho += rho_step)
		for (double ct = ct_from; ct <= ct_to; ct += ct_step)
			thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call(to_compute);
	std::vector<std::vector<double>> data = thead_core->get_value();

	Superposition medium_type = this->global_conf->medium_superposition();
	std::vector<std::vector<double>> plot_data;
	for (auto i : data) {
		double magnitude = 0;
		for (std::size_t j = 4; j < i.size(); j++)
			magnitude += i[j];
		plot_data.push_back({magnitude,i[0],i[1]});
	}

	GnuPlot* plot = new GnuPlot(this->global_conf->gnp_script_path());
	plot->set_gnuplot_bin(this->global_conf->path_gnuplot_binary());
	plot->set_colormap(this->global_conf->plot_color_map());
	plot->set_ox_label("ct, m");
	plot->set_oy_label("rho, m");
	plot->set_oz_label("Ex, V/m");
	plot->grid_on();
	plot->cage_on();
	plot->plot3d(plot_data);
	if ( this->global_conf->call_gnuplot() ) plot->call_gnuplot();
}

void PlotModel::__from_ct_z (const std::vector<std::pair<Component,AbstractField*>>& to_compute)
{
	std::size_t thread_num = this->global_conf->thread_number();
	double ct_from = this->global_conf->receiver_vt()[0];
	double ct_step = this->global_conf->receiver_vt()[1];
	double ct_to = this->global_conf->receiver_vt()[2];
	double rho = this->global_conf->receiver_rho()[0];
	double phi = this->global_conf->receiver_phi()[0];
	double z_from = this->global_conf->receiver_z()[0];
	double z_step = this->global_conf->receiver_z()[1];
	double z_to = this->global_conf->receiver_z()[2];

	Manager* thead_core;
	if (!this->global_conf->safe_mode()) thead_core = new Manager( thread_num );
	else thead_core = new SafeManager( thread_num, this->global_conf );
	thead_core->progress_bar( this->global_conf->print_progress() );

	for (double z = z_from; z <= z_to; z += z_step)
		for (double ct = ct_from; ct <= ct_to; ct += ct_step)
			thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call(to_compute);
	std::vector<std::vector<double>> data = thead_core->get_value();

	Superposition medium_type = this->global_conf->medium_superposition();
	std::vector<std::vector<double>> plot_data;
	for (auto i : data) {
		double magnitude = 0;
		for (std::size_t j = 4; j < i.size(); j++)
			magnitude += i[j];
		plot_data.push_back({magnitude,i[0],i[3]});
	}

	GnuPlot* plot = new GnuPlot( this->global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( this->global_conf->path_gnuplot_binary() );
	plot->set_colormap(this->global_conf->plot_color_map());
	plot->set_ox_label("ct, m");
	plot->set_oy_label("z, m");
	plot->set_oz_label("Ex, V/m");
	plot->grid_on();
	plot->cage_on();
	plot->plot3d(plot_data);
	if ( this->global_conf->call_gnuplot() ) plot->call_gnuplot();
}

void PlotModel::__from_x_y (const std::vector<std::pair<Component,AbstractField*>>& to_compute)
{
	throw std::logic_error("Not implemented!");
}

/* void ReadyModel::__from_rho_phi (double ct, double z)
{
	// not tested and not used
	float R = Config::plane_disk_radius();
	float A0 = Config::plane_disk_magnitude();
	float eps_r = Config::plane_disk_epsr();
	float mu_r = Config::plane_disk_mur();

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	linear->set_yterms_num( Config::magnetic_term_num() );

	std::size_t thread_num = Config::thread_number();
	Manager* thead_core = new Manager( thread_num );
	thead_core->progress_bar( Config::print_progress() );

	for (double rho = 0; rho <= 2; rho += 0.01)
		for (double phi = 0; phi <= 2 * M_PI; phi += 0.01)
			thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call( *linear, "electric_x" );
	std::stack<std::vector<double>> res = thead_core->get_value();

	std::vector<std::tuple<double, double, double>> plot_data;
	while (!res.empty()) {
		std::vector<double> tmp = res.top(); res.pop();
		auto point = std::make_tuple(tmp[2], tmp[3], tmp[0]);
		plot_data.push_back(point);
	}


	GnuPlot* plot = new GnuPlot( Config::gnp_script_path() );
	plot->set_gnuplot_bin( Config::path_gnuplot_binary() );
	plot->set_ox_label("rho, m");
	plot->set_oy_label("phi, m");
	plot->set_oz_label("Ex, V/m");
	plot->grid_on();
	plot->cage_on();
	plot->plot3d(plot_data);
	if ( Config::call_gnuplot() ) plot->call_gnuplot();
} */

/* void ReadyModel::missile_effect_length ()
{
	std::vector<std::vector<double>> plot_data;
	for (double z = 1; z <= 20; z += 0.1) {
		double val1 = std::sqrt(2*2 + z*z) - z; // red
		double val2 = std::sqrt(5*5 + z*z) - z; // blue
		double val3 = 2*2 / z; // black
		double val4 = 2*5 / z; // green
		auto point = {z, val1, val2, val3, val4};
		plot_data.push_back(point);
	}

	GnuPlot* plot = new GnuPlot( Config::gnp_script_path() );
	plot->set_gnuplot_bin( Config::path_gnuplot_binary() );
	plot->set_ox_label("z, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	// plot->set_logscale_ox(true);
	std::vector<std::string> title = {"c * tau(R=0.7), m", 
		"c * tau(R=2), m", "0.7*0.7 / (2z), m", "2*2 / z, m"};
	plot->plot_multi(plot_data, title);
	if ( Config::call_gnuplot() ) plot->call_gnuplot();
} */

/* void ReadyModel::magnetic_static_magnitude (double eps, double rho, double phi)
{
	// not tested and not used
	float R = Config::plane_disk_radius();
	float A0 = Config::plane_disk_magnitude();
	float eps_r = Config::plane_disk_epsr();
	float mu_r = Config::plane_disk_mur();

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	linear->set_yterms_num(70);

	std::vector<std::pair<double,double>> plot_data;
	for (double z = 0; z <= 5; z += 0.1) {
		double As = linear->static_magnitude(rho,phi,z,eps);
		// double As = linear->static_magnitude(z);
		auto point = std::make_pair(z, As);
		plot_data.push_back(point);
	}

	GnuPlot* plot = new GnuPlot( Config::gnp_script_path() );
	plot->set_gnuplot_bin( Config::path_gnuplot_binary() );
	plot->set_ox_label("z, m");
	plot->set_oy_label("As (Hy), A/m");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(plot_data);
	if ( Config::call_gnuplot() ) plot->call_gnuplot();
} */

/* void ReadyModel::magnetic_static_terms (double rho, double phi, double z)
{
	float R = Config::plane_disk_radius();
	float A0 = Config::plane_disk_magnitude();
	float eps_r = Config::plane_disk_epsr();
	float mu_r = Config::plane_disk_mur();

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);

	MissileField* linear0 = new MissileField(source, medium);
	linear0->set_yterms_num(0);
	auto term0 = [&linear0, rho, phi, z] (double arg) {
		return linear0->magnetic_y(arg, rho, phi, z);
	};

	MissileField* linear1 = new MissileField(source, medium);
	linear1->set_yterms_num(1);
	auto term1 = [&linear1, rho, phi, z] (double arg) {
		return linear1->magnetic_y(arg, rho, phi, z);
	};

	MissileField* linear5 = new MissileField(source, medium);
	linear5->set_yterms_num(5);
	auto term5 = [&linear5, rho, phi, z] (double arg) {
		return linear5->magnetic_y(arg, rho, phi, z);
	};

	MissileField* linear20 = new MissileField(source, medium);
	linear20->set_yterms_num(20);
	auto term20 = [&linear20, rho, phi, z] (double arg) {
		return linear20->magnetic_y(arg, rho, phi, z);
	};

	MissileField* linear40 = new MissileField(source, medium);
	linear40->set_yterms_num(40);
	auto term40 = [&linear40, rho, phi, z] (double arg) {
		return linear40->magnetic_y(arg, rho, phi, z);
	};

	MissileField* linear100 = new MissileField(source, medium);
	linear100->set_yterms_num(100);
	auto term100 = [&linear100, rho, phi, z] (double arg) {
		return linear100->magnetic_y(arg, rho, phi, z);
	};

	std::size_t thread_num = Config::thread_number();
	Manager* thead_core = new Manager( thread_num );
	thead_core->progress_bar( Config::print_progress() );

	for (double ct = 0.1; ct <= 100; ct += 0.01)
		thead_core->add_argument( {ct} );

	thead_core->call( {term0, term1, term5, term20, term40, term100} );
	std::stack<std::vector<double>> res = thead_core->get_value();

	std::vector<std::vector<double>> plot_data;
	while (!res.empty()) {
		std::vector<double> tmp = res.top(); res.pop();
		plot_data.push_back(tmp);
	}

	typedef std::vector<double> vector;
	std::sort(plot_data.begin(), plot_data.end(),
	[] (const vector &a, const vector &b) -> bool
	{ 
		return a.front() < b.front();
	});

	GnuPlot* plot = new GnuPlot( Config::gnp_script_path() );
	plot->set_gnuplot_bin( Config::path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	plot->set_logscale_ox(true);
	std::vector<std::string> title = { "N = 0", 
									   "N = 1", 
									   "N = 5", 
									   "N = 20", 
									   "N = 40", 
									   "N = 100"
	};
	plot->plot_multi(plot_data, title);
	if ( Config::call_gnuplot() ) plot->call_gnuplot();
} */
