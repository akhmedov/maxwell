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
	this->global_conf = conf;
	this->model_pointer = {
		&PlotModel::__Ex_from_ct, 
		&PlotModel::__Hy_from_ct
	};
}

void PlotModel::call (const PlotModel::Name& model_name)
{
	this->model_pointer[model_name - 1](*this);
}

void PlotModel::__Ex_from_ct ()
{
	float R = this->global_conf->plane_disk_radius();
	float A0 = this->global_conf->plane_disk_magnitude();
	float eps_r = this->global_conf->plane_disk_epsr();
	float mu_r = this->global_conf->plane_disk_mur();

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	linear->set_yterms_num( this->global_conf->magnetic_term_num() );

	std::size_t thread_num = this->global_conf->thread_number();
	Manager* thead_core = new Manager( thread_num );
	thead_core->progress_bar( this->global_conf->print_progress() ); 

	double ct_from = this->global_conf->receiver_vt()[0];
	double ct_step = this->global_conf->receiver_vt()[1];
	double ct_to = this->global_conf->receiver_vt()[2];
	double rho = this->global_conf->receiver_rho()[0];
	double phi = this->global_conf->receiver_phi()[0];
	double z = this->global_conf->receiver_z()[0];

	for (double ct = ct_from; ct <= ct_to; ct += ct_step)
		thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call( *linear, "electric_x" );
	std::stack<std::vector<double>> res = thead_core->get_value();

	std::vector<std::pair<double, double>> plot_data;
	while (!res.empty()) {
		std::vector<double> tmp = res.top(); res.pop();
		auto point = std::make_pair(tmp[1], tmp[0]);
		plot_data.push_back(point);
	}

	// sort of output
	// TODO: use std::set for online sorting
	typedef std::pair<double, double> pair;
	std::sort (plot_data.begin(), plot_data.end(), 
		[] (const pair& i, const pair& j) -> bool 
		{ return i.first < j.first; }
	);

	GnuPlot* plot = new GnuPlot( this->global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( this->global_conf->path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("Ex, V/m");
	plot->grid_on();
	plot->cage_on();
	plot->set_logscale_ox(true);
	plot->plot2d(plot_data);
	if ( this->global_conf->call_gnuplot() ) plot->call_gnuplot();
}

void PlotModel::__Hy_from_ct ()
{
	float R = this->global_conf->plane_disk_radius();
	float A0 = this->global_conf->plane_disk_magnitude();
	float eps_r = this->global_conf->plane_disk_epsr();
	float mu_r = this->global_conf->plane_disk_mur();

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	linear->set_yterms_num( this->global_conf->magnetic_term_num() );

	std::size_t thread_num = this->global_conf->thread_number();
	Manager* thead_core = new Manager( thread_num );
	thead_core->progress_bar( this->global_conf->print_progress() ); 

	double ct_from = this->global_conf->receiver_vt()[0];
	double ct_step = this->global_conf->receiver_vt()[1];
	double ct_to = this->global_conf->receiver_vt()[2];
	double rho = this->global_conf->receiver_rho()[0];
	double phi = this->global_conf->receiver_phi()[0];
	double z = this->global_conf->receiver_z()[0];

	for (double ct = ct_from; ct <= ct_to; ct += ct_step)
		thead_core->add_argument( {ct,rho,phi,z} );

	thead_core->call( *linear, "magnetic_y" );
	std::stack<std::vector<double>> res = thead_core->get_value();

	std::vector<std::pair<double, double>> plot_data;
	while (!res.empty()) {
		std::vector<double> tmp = res.top(); res.pop();
		auto point = std::make_pair(tmp[1], tmp[0]);
		plot_data.push_back(point);
	}

	typedef std::pair<double, double> pair;
	std::sort (plot_data.begin(), plot_data.end(), 
		[] (const pair& i, const pair& j) -> bool 
		{ return i.first < j.first; }
	);

	GnuPlot* plot = new GnuPlot( this->global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( this->global_conf->path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("Hy, A/m");
	plot->grid_on();
	plot->cage_on();
	plot->set_logscale_ox(true);
	plot->plot2d(plot_data);
	if ( this->global_conf->call_gnuplot() ) plot->call_gnuplot();
}
