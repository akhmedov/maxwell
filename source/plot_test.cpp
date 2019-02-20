//
//  plot_test.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 31.01.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "plot_test.hpp"

Config* PlotTest::global_conf = new Config();

int main()
{
	PlotTest::set_options();

	/* DERIVATIVE OF I1 AND I2 */

	/* std::cout << std::endl << "PlotTest::";
	std::cout << "I1_time_partder(rho = 0.5, z = 0.6, R = 1) " << std::endl;
	PlotTest::I1_time_partder (0.5,0.6);
	
	std::cout << std::endl << "PlotTest::";
	std::cout << "I1_time_partder(rho = 1, z = 0.6, R = 1) " << std::endl;
	PlotTest::I1_time_partder (1,0.6);
	
	std::cout << std::endl << "PlotTest::";
	std::cout << "I1_time_partder(rho = 2, z = 1, R = 1) " << std::endl;
	PlotTest::I1_time_partder (2,1);

	std::cout << std::endl << "PlotTest::";
	std::cout << "I2_time_partder(rho = 0.5, z = 0.6, R = 1) " << std::endl;
	PlotTest::I2_time_partder (0.5,0.6);
	
	std::cout << std::endl << "PlotTest::";
	std::cout << "I2_time_partder(rho = 1, z = 0.6, R = 1) " << std::endl;
	PlotTest::I2_time_partder (1,0.6);
	
	std::cout << std::endl << "PlotTest::";
	std::cout << "I2_time_partder(rho = 2, z = 1, R = 1) " << std::endl;
	PlotTest::I2_time_partder (2,1);

	std::cout << std::endl << "PlotTest::";
	std::cout << "I2_time_partder(rho = 0.5, z = 1.21, R = 0.1) " << std::endl;
	PlotTest::I2_time_partder (2,1,0.1); */

	/* I1 NUMERICAL INTEGRAL WITH N POINTS */

	/* std::cout << std::endl << "PlotTest::";
	std::cout << "I1_numeric_integral(1e3,1e2) " << std::endl;
	PlotTest::I1_numeric_integral(1e3,1e2);

	std::cout << std::endl << "PlotTest::";
	std::cout << "I1_numeric_integral(1e3,1e2) " << std::endl;
	PlotTest::I1_numeric_integral(1e3,1e1);

	std::cout << std::endl << "PlotTest::";
	std::cout << "I1_numeric_integral(1e3,1e2) " << std::endl;
	PlotTest::I1_numeric_integral(1e2,1e1);

	std::cout << std::endl << "PlotTest::";
	std::cout << "I1_numeric_integral(1e3,1e2) " << std::endl;
	PlotTest::I1_numeric_integral(2e1,2e1); */

	/* NONLINEAR ELECTRIC CURRENT DISTRIBUTION J' */

	/* std::cout << std::endl << "PlotTest::";
	std::cout << "nonlinear_current(phi=0,z=0) " << std::endl;
	PlotTest::nonlinear_current(3.14159,0.5); */

	/* KERR TE WAVE UNDERINTEGRAL EXPRESION FROM NU */

	// std::cout << std::endl << "PlotTest::";
	// std::cout << "kerr_under_integral_from_nu(r,r',R=1) " << std::endl;
	// std::vector<double> r = {3, 0, 0, std::sqrt(2)};
	// std::vector<double> r_perp = {std::sqrt(3), 0.5, 0, std::sqrt(1.5)};
	// PlotTest::kerr_under_integral_from_nu(r, r_perp);

	/* DUHAMEL LINEAR FIELD */

	// std::cout << std::endl << "PlotTest::";
	// std::cout << "plot_source() " << std::endl;
	// PlotTest::plot_source(1);
	// PlotTest::gauss_signal(1);

	// std::cout << std::endl << "PlotTest::";
	// std::cout << "Ex_derivative() " << std::endl;
	// PlotTest::Ex_derivative(0.1);

	/* ENERGY PLOT */

	std::cout << std::endl << "PlotTest::";
	std::cout << "plot_energy_slyse() " << std::endl;
	// PlotTest::plot_energy_slyse(1, 0.50);
	// PlotTest::plot_energy_slyse(1, 1.00);
	// PlotTest::plot_energy_slyse(1, 2.50);
	// PlotTest::plot_energy_slyse(1, 5.00);
	// PlotTest::plot_energy_slyse(1,10.00);
	// PlotTest::plot_energy_slyse(1,15.00);
	// PlotTest::plot_energy_slyse(1,20.00);
	// PlotTest::plot_energy_slyse(1,30.00);
	// PlotTest::plot_energy_slyse(1,40.00);

	PlotTest::plot_energy_distribution(1, 5);

	/* std::cout << std::endl << "PlotTest::plot_energy_max(tau) ... " << std::endl;
	PlotTest::plot_energy_max();
	std::cout << "Done." << std::endl; */

	/* PLOT FIELD FROM VT FOR TAU1 TAU2 TAU3 */

	// std::cout << "PlotTest::plot_Ex_from_tau" << std::endl;
	// PlotTest::plot_Ex_from_tau();
	// PlotTest::recive_emp(0, 0, 2);

	// PlotTest::energy_compare(1,2);
	// PlotTest::energy_iterfer_sinc();

	/* NOISE POWER */

	// PlotTest::awgn_power({100,1000});
}

void PlotTest::set_options ()
{
	PlotTest::global_conf->path_gnuplot_binary("gnuplot/bin/gnuplot");
	/* TODO: BUG - not used config param
	PlotTest::global_conf->plot_color_map(Colormap::parula); */
}

void PlotTest::gauss_signal (std::size_t order)
{
	double ct1 = 1.7, ct2 = 4, rho = 0, phi = 0, z = 2;
	double R = 1, A0 = 1, tau = 1, eps_r = 1, mu_r = 1;

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau,order] (double vt) {return Function::gauss_perp(vt,tau,order);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	for (double vt = ct1; vt < ct2; vt += 0.001) {
		line = {vt, duhamel->electric_x(vt, rho, phi, z)};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("f(ct)");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(plot_data);
	plot->call_gnuplot();
}

void PlotTest::awgn_power (const std::vector<std::size_t>& samples)
{
	std::vector<std::vector<double>> data;
	std::vector<double> line;

	for  (double sigma = 0.1; sigma < 10; sigma += 0.1) {
		line = {sigma};
		AdditiveWhiteGaussian* noise = new AdditiveWhiteGaussian(0,sigma);
		for (auto i : samples) {
			double Pn = noise->power(0,0,0,i);
			line.push_back(Pn);
		}
		data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( "power.gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("sigma");
	plot->set_oy_label("Pn, V*V");
	plot->grid_on(false);
	plot->cage_on();
	std::vector<std::string> title;
	for (auto i : samples) title.push_back("N = " + std::to_string(i));
	plot->plot_multi(data, title);
}

void PlotTest::emp_duration (double rho, double tau0)
{
	std::vector<std::vector<double>> data;
	std::vector<double> line;

	for (double z = 0; z < 10; z += 0.01) {
		double tau = serial::dataset::updisk_emp_duraton(tau0,rho,z);
		line = {z, tau};
		data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("z, m");
	plot->set_oy_label("tau, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(data);
	plot->call_gnuplot();
}

void PlotTest::energy_iterfer_sinc()
{
	double R = 1, A0 = 1, tau = 0.5;
	double eps_r = 1, mu_r = 1;
	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);
	auto property = &AbstractField::energy_cart;
	auto compute = std::make_pair(property, (AbstractField*)duhamel);

	Manager<2>* thead_core = new Manager<2>(4, NULL);
	thead_core->progress_bar(true);
	for (double z = 8; z < 16; z += 0.1) {
		double rho = 0;
		double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
		double to = tau + std::sqrt((rho+R)*(rho+R) + z*z);
		thead_core->add_argument( {0,0,z,from,to} );
	}

	thead_core->call({compute});
	std::vector<std::vector<double>> data = thead_core->get_value();

	for (auto&& i : data) {
		i.erase(i.begin(), i.begin() + 2); // erase x y
		// i[1] *= i[0] * i[0];
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("z, m");
	plot->set_oy_label("W, V*V/m2");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(data);
	plot->call_gnuplot();
}

void PlotTest::energy_compare (double tau1, double tau2)
{
	double R = 1, A0 = 1;
	auto field = [R,A0] (double tau) {
		double eps_r = 1, mu_r = 1;
		Homogeneous* medium = new Homogeneous(mu_r, eps_r);
		UniformPlainDisk* source = new UniformPlainDisk(R, A0);
		MissileField* linear = new MissileField(source, medium);
		FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
		free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
		LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);
		auto property = &AbstractField::energy_cart;
		return std::make_pair(property, (AbstractField*)duhamel);
	};

	for (double x = 0; x < 3; x += 0.1) {
		double z = 0;
		double tau = std::max(tau1,tau2);
		double from = (x > R) ? std::sqrt((x-R)*(x-R) + z*z) : z;
		double to = tau + std::sqrt((x+R)*(x+R) + z*z);
		Manager<2>* thead_core = new Manager<2>(1, NULL);
		thead_core->progress_bar(false);
		thead_core->add_argument( {x,0,z,from,to} );
		thead_core->call({field(tau1), field(tau2)});
		std::vector<std::vector<double>> data = thead_core->get_value();
		
		std::cout << data[0][0] << " : "  << data[0][3] << " , " << data[0][4] << std::endl;
	}
}

void PlotTest::recive_emp (double rho, double phi, double z)
{
	double R = 1, A0 = 1, tau = 0.5;
	double eps_r = 1, mu_r = 1;

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::sinc(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	auto emp = [duhamel, rho, phi, z] (double vt) {
		return duhamel->electric_x(vt, rho, phi, z);
	};

	for (double vt = z + tau - 0.5; vt < z + 2 * tau; vt += 0.0001) {
		// line = {vt, emp(vt)};
		line = {vt, Math::derivat4(emp, vt)};
		std::cout << line[0] << ' ' << line[1] << std::endl;
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("f(ct)");
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(plot_data);
	plot->call_gnuplot();
}

void PlotTest::plot_log_from_n ()
{
	std::vector<std::vector<double>> data;
	for (int n = 1; n <= 16; n++) {
		std::vector<double> line = {(float)n};
		line.push_back(Math::log(n,1000));
		data.push_back(line);
	}

	GnuPlot* plot = new GnuPlot( "log_n.gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("n");
	plot->set_oy_label("log_n(1e3)");
	plot->grid_on(false);
	plot->plot2d(data);
}

void PlotTest::plot_Ex_from_tau ()
{
	auto field = [] (double tau) {
		double R = 1, A0 = 1;
		double eps_r = 1, mu_r = 1;
		Homogeneous* medium = new Homogeneous(mu_r, eps_r);
		UniformPlainDisk* source = new UniformPlainDisk(R, A0);
		MissileField* on = new MissileField(source, medium);
		auto property = &AbstractField::electric_x;
		if (tau == 0) {
			MeandrPeriod* rect = new MeandrPeriod(R, A0, 0.7);
			SquaredPulse* meandr = new SquaredPulse(rect, medium);
			return std::make_pair(property, (AbstractField*)meandr);
		}
		FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
		free_shape->set_time_depth([tau] (double vt) {return Function::smoozed_rect(vt,0.5,tau);});
		LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, on, NULL);
		return std::make_pair(property, (AbstractField*)duhamel);
	};

	Manager<0>* thead_core = new Manager<0>(4, NULL);
	thead_core->progress_bar(true);

	for (double vt = 1.8; vt <= 3; vt += 0.01)
			thead_core->add_argument( {vt,0,0,2} );

	thead_core->call({field(0), field(0.1), field(0.2)});
	std::vector<std::vector<double>> data = thead_core->get_value();

	for (auto&& i : data)
		i.erase(i.begin() + 1, i.begin() + 4); // erase rho phi z

	GnuPlot* plot = new GnuPlot( "Ex_from_tau.gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("vt, m");
	plot->set_oy_label("Ex, V/m");
	// plot->grid_on(false);
	// plot->cage_on();
	std::vector<std::string> title = {"tau = 0 vt",
									  "tau = 0.1 vt", 
									  "tau = 0.2 vt"};
	plot->plot_multi(data, title);
}

void PlotTest::plot_energy_max ()
{

	PlotTest::global_conf->field_component(FieldComponent::W);
	PlotTest::global_conf->impulse_shape(ImpulseShape::gauss);
	// PlotTest::global_conf->duration(0);
	
	double R = 1, A0 = 1;

	auto field = [R,A0] (double tau) {
		double eps_r = 1, mu_r = 1;
		Homogeneous* medium = new Homogeneous(mu_r, eps_r);
		UniformPlainDisk* source = new UniformPlainDisk(R, A0);
		MissileField* linear = new MissileField(source, medium);
		FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
		free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
		LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);
		auto property = &AbstractField::energy_cart;
		return std::make_pair(property, duhamel);
	};

	Manager<2>* thead_core = new Manager<2>(4, NULL);
	thead_core->progress_bar(true);

	for (double z = 0.5; z <= 15; z += 0.1) {
		double rho = 0;
		double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
		double to = 1 + std::sqrt((rho+R)*(rho+R) + z*z);
		thead_core->add_argument( {0,0,z,from,to} );
	}

	thead_core->call({field(0.5), field(1)});
	std::vector<std::vector<double>> data = thead_core->get_value();
	
	double max0 = 0;
	for (auto i : data) 
		if (std::abs(i[2] - 15) < 1e-5) 
			max0 = i[3] * i[2] * i[2];

	for (auto&& i : data) {
		i[3] *= i[2] * i[2] / max0; // norm W
		i[4] *= i[2] * i[2] / max0; // norm W
		i.erase(i.begin(), i.begin() + 2); // erase x and y
	}

	GnuPlot* plot = new GnuPlot( "onaxis_w_z2.gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("z, m");
	plot->set_oy_label("W, V*V");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"tau = 0.5 vt", 
									  "tau = 1.0 vt"};
	plot->plot_multi(data, title);
}

void PlotTest::plot_energy_distribution (double tau, double max_z)
{
	auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

	double R = 1, A0 = 1;
	double eps_r = 1, mu_r = 1;

	double x = 0;

	PlotTest::global_conf->field_component(FieldComponent::W);
	PlotTest::global_conf->impulse_shape(ImpulseShape::gauss);
	PlotTest::global_conf->duration(tau);

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::sinc(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

	Manager<0>* thead_core = new Manager<0>(6, NULL);
	thead_core->progress_bar(true);

	for (double y = -2*R; y <= 2*R; y += 0.05) {
		for (double z = 0; z <= max_z; z += 0.05) {
			double rho = std::sqrt(x*x + y*y);
			double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau + std::sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to+0.01} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = std::make_pair(property, linear);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	GnuPlot* plot = new GnuPlot( str_of(tau) + "_" + str_of(max_z) + ".gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("y, m");
	plot->set_oy_label("z, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 1, 2);
}

void PlotTest::plot_energy_slyse (double tau, double z)
{
	auto str_of = [] (double a) {return std::to_string(a).substr(0,5);};

	double R = 1, A0 = 1;
	double eps_r = 1, mu_r = 1;
	double range = z/2;

	PlotTest::global_conf->field_component(FieldComponent::W);
	PlotTest::global_conf->impulse_shape(ImpulseShape::sinc);
	PlotTest::global_conf->duration(tau);

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source);
	free_shape->set_time_depth([tau] (double vt) {return Function::sinc(vt,tau);});
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, linear, NULL);

	Manager<0>* thead_core = new Manager<0>(6, NULL);
	thead_core->progress_bar(true);

	for (double x = -range; x <= range; x += 0.05) {
		for (double y = -range; y <= range; y += 0.05) {
			double rho = std::sqrt(x*x + y*y);
			double from = (rho > R) ? std::sqrt((rho-R)*(rho-R) + z*z) : z;
			if (from - 0.01 > 0) from -= 0.01;
			double to = tau + std::sqrt((rho+R)*(rho+R) + z*z);
			thead_core->add_argument( {x,y,z,from,to+0.01} );
		}
	}

	auto property = &AbstractField::energy_cart;
	auto function = std::make_pair(property, linear);
	// auto function = std::make_pair(property, duhamel);
	thead_core->call({function});

	std::vector<std::vector<double>> data = thead_core->get_value();

	std::vector<double> max0 = data[0];
	for (auto point : data)
		if (point[3] > max0[3])
			max0 = point;

	// for (auto&& i : data) {
	// 	i[3] *= z*z / max0[3]; // norm W
	// 	i.erase(i.begin()+2,i.begin()+5); // erase z, from, to
	// }

	std::cout << "Wmax (" << str_of(max0[0]) << ',' << str_of(max0[1]) << ',' << str_of(max0[2]) << ") = " << str_of(max0[3]) << std::endl;

	/* std::vector<std::vector<double>> agumented;
	for (auto&& i : data) {
		double x = i[0], y = i[1], W = i[2];
		if (x > 1e-8) agumented.push_back( {-x,y,W} );
		if (y > 1e-8) agumented.push_back( {x,-y,W} );
		if (x > 1e-8 && y > 1e-8) agumented.push_back( {-x,-y,W} );
	}
	data.insert(std::end(data), std::begin(agumented), std::end(agumentat)); */

	GnuPlot* plot = new GnuPlot( str_of(tau) + "_" + str_of(z) + ".gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("x, m");
	plot->set_oy_label("y, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(data, 0, 1);
}

void PlotTest::Ex_derivative (double tau)
{
	double R = 1;
	double eps_r = 1;
	double mu_r = 1;
	double xi3 = 0;
	double A0 = 1;

	Homogeneous* linear_medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, linear_medium);

	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	auto f = [linear, tau] (double vt) {
		return linear->electric_x(vt - tau, 1, 0, 2);
	};

	for (double vt = 1.7; vt < 3; vt += 0.01) {
		double anal = 0;
		double numerics = Math::derivat3(f, vt);

		line = {vt, anal, numerics};
		plot_data.push_back(line);
		line.clear();
	}


}

void PlotTest::plot_source (std::size_t order)
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	// auto f = [] (double x) { return Function::rect(x, 1); };
	// auto f = [] (double x) { return Function::sin(x, 1); };
	// auto f = [] (double x) { return Function::sinc(x, 1); };
	// auto f = [] (double x) { return Function::gauss(x, 1); };
	// auto f = [] (double x) { return Function::smoozed_rect(x, 1, 0.2); };
	// auto f = [] (double x) { return Function::gauss(x, 1); };
	auto f = [order] (double x) { return Function::gauss_perp_normed(x, 1, order); };
	// auto f = [] (double x) { return Function::gauss(x, 2) * Function::sinc(x, 2, 10); };

	for (double vt = -1; vt < 3; vt += 0.01) {
		double source = f(vt);
		line = {vt, source};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("ct, m");
	plot->set_oy_label("f(t)"+std::to_string(order));
	plot->grid_on();
	plot->cage_on();
	plot->plot2d(plot_data);
	plot->call_gnuplot();
}


void PlotTest::kerr_under_integral_from_nu (const std::vector<double> &r, const std::vector<double> &r_perp, double R)
{
	double vt = r[0];
	double vt_perp = r_perp[0];
	double rho = r[1];
	double rho_perp = r_perp[1];
	double phi = r[2];
	double phi_perp = r_perp[2];  UNUSED(phi_perp);
	double z = r[3];
	double z_perp = r_perp[3];

	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	double max_vt = vt - z + z_perp; // grater then zero
	double delta_vt = vt - vt_perp;
	double delta_z = z - z_perp;
	double casual = std::sqrt(delta_vt*delta_vt - delta_z*delta_z);
	if (std::isnan(casual)) throw std::invalid_argument("Casuality princip is not supported at r and r'");
	double max_nu =  PERIODS_NU * std::abs(rho - rho_perp + casual);
	if (max_nu <= 1e-5) throw std::invalid_argument("max_nu is zero value for current r and r'");

	for (double nu = 0; nu < 3 * max_nu; nu += 0.001) {

		double term1 = 0;
		for (int m = -3; m <= 3; m += 2)
			term1 += KerrAmendment::x_trans( m, nu, rho, phi) *
				     KerrAmendment::N_sum(R, m, nu, max_vt, rho_perp, z_perp);

		double bessel = jn(0,nu * casual) + jn(2,nu * casual);
		double term2 = 0;
		for (int m = -3; m <= 3; m += 2)
			term2 += KerrAmendment::x_trans( m, nu, rho, phi) *
				     KerrAmendment::N_sum(R, m, nu, vt_perp, rho_perp, z_perp);
		term2 *= nu * nu * delta_vt * bessel / 2;

		line = {nu, term1, term2, term1 + term2};
		plot_data.push_back(line);
		line.clear();

		std::cout << nu << ' ' << term1 << ' ' << term2 << std::endl;
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	// plot->set_colormap(Colormap::parula);
	plot->set_ox_label("nu");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"delta term",
									  "step term",
									  "sum"};
	plot->plot_multi(plot_data, title);
	plot->call_gnuplot();
}

void PlotTest::nonlinear_current (double phi, double z)
{
	double R = 1;
	double eps_r = 1;
	double mu_r = 1;
	double xi3 = 0;
	// double inv_em_relation = NonlinearMedium::EPS0 * eps_r / NonlinearMedium::MU0 * mu_r;
	double A0 = 2;//* inv_em_relation;

	Homogeneous* linear_medium = new Homogeneous(mu_r, eps_r);
	KerrMedium* kerr_medium = new KerrMedium(mu_r, eps_r, xi3, 0);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, linear_medium);
	KerrAmendment* non_linear = new KerrAmendment(linear, kerr_medium, source);
	std::vector<std::vector<double>> plot_data;

	for (double rho = 0.0; rho <= 1.5; rho += 0.01) {
		for (double ct = 0.4; ct <= 1.5; ct += 0.01) {
			if (ct - z <= 0) {
				plot_data.push_back({0,ct,rho});
			} else {
				double j1 = non_linear->current_x(ct,rho,phi,z);
				plot_data.push_back({j1,ct,rho});
			}
		}
	}

	GnuPlot* plot = new GnuPlot(PlotTest::global_conf->gnp_script_path());
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	// plot->set_colormap(Colormap::parula);
	plot->set_ox_label("ct, m");
	plot->set_oy_label("rho, m");
	plot->set_oz_label("");
	plot->grid_on();
	plot->cage_on();
	plot->plot3d(plot_data);
	plot->call_gnuplot();
}

void PlotTest::I1_numeric_integral (std::size_t points, double limit)
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	double a = 1, b = 0.5;

	for (double c = 0; c < 2; c += 0.05) {

		auto f = [a, b, c] (double x) {
			if (b == 0) return a * 0.5 * j1(a*x) * j0(c*x);
			if (x == 0) return 0.0;
			return a * j1(a*x) * j1(b*x) * j0(c*x) / (x * b);
		};

		double anal = MissileField::int_bessel_011(c, b, a);
		Simpson I = Simpson(points);
		double numerics = I.value(0, limit, f);

		line = {c, anal, numerics};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::parula);
	plot->set_ox_label("ct, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"analitical integral",
									  "numerical aproximathin"};
	plot->plot_multi(plot_data, title);
	plot->call_gnuplot();
}

void PlotTest::numeric_perp(double omega)
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	for (double x = 0; x <= 2 * M_PI; x += 0.01) {

		auto sin = [omega] (double x) {
			return std::sin(omega * x);
		};
		
		double anal = std::cos(omega * x);
		double num3 = Math::derivat3(sin, omega * x);
		double num4 = Math::derivat4(sin, omega * x);

		line = {omega*x, sin(omega*x), anal, num3, num4};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::parula);
	plot->set_ox_label("ct, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"test function",
									  "analitical derivative",
									  "3th order aproximathin",
									  "4th order aproximathin"};
	plot->plot_multi(plot_data, title);
	plot->call_gnuplot();
}

void PlotTest::I1_I2_versus (double rho, double z, double R)
{
	UNUSED(rho); UNUSED(z); UNUSED(R);
	throw std::logic_error("PlotTest::I1_I2_versus() is not implemented");
}

void PlotTest::I1_time_partder (double rho, double z, double R)
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	for (double ct = z+0.01; ct <= z + 3; ct += 0.01) {

		auto I1 = [rho, R, z] (double ct) {
			double vt_z = std::sqrt(ct * ct - z * z);
			return MissileField::int_bessel_011(vt_z,rho,R);
		};
		
		double anal = KerrAmendment::int_bessel_011_perp(ct,z,rho,R);
		double num3 = Math::derivat3(I1, ct);
		double num4 = Math::derivat4(I1, ct);

		line = {ct, I1(ct), anal, num3, num4};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::gray);
	plot->set_ox_label("ct, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"I1 form time dependency",
									  "alalitical derivative", 
									  "3th order aproximathin",
									  "4th order aproximathin"};
	plot->plot_multi(plot_data, title);
	plot->call_gnuplot();
}

void PlotTest::I2_time_partder (double rho, double z, double R)
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	for (double ct = z+0.1; ct <= z + 3; ct += 0.01) {

		auto I2 = [rho, R, z] (double ct) {
			double vt_z = std::sqrt(ct * ct - z * z);
			return MissileField::int_bessel_001(vt_z,rho,R);
		};
		
		double anal = KerrAmendment::int_bessel_001_perp(ct,z,rho,R);
		double num3 = Math::derivat3(I2, ct);
		double num4 = Math::derivat4(I2, ct);

		line = {ct, I2(ct), anal, num3, num4};
		plot_data.push_back(line);
		line.clear();
	}

	GnuPlot* plot = new GnuPlot( PlotTest::global_conf->gnp_script_path() );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_colormap(Colormap::parula);
	plot->set_ox_label("ct, m");
	plot->set_oy_label("");
	plot->grid_on();
	plot->cage_on();
	std::vector<std::string> title = {"I2 form time dependency",
									  "alalitical derivative", 
									  "3th order aproximathin",
									  "4th order aproximathin"};
	plot->plot_multi(plot_data, title);
	plot->call_gnuplot();
}

void PlotTest::kerr_ammend_undeintegral (double ct_perp, 
										 double phi_perp, 
									     double rho_perp, 
										 double z_perp, 
										 double R)
{
	UNUSED(ct_perp); UNUSED(phi_perp); UNUSED(rho_perp); UNUSED(z_perp); UNUSED(R);
	throw std::logic_error("Not implemented!");
}
