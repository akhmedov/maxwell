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
	// PlotTest::plot_source();

	// std::cout << std::endl << "PlotTest::";
	// std::cout << "Ex_derivative() " << std::endl;
	// PlotTest::Ex_derivative(0.1);

	/* ENERGY PLOT */

	// std::cout << std::endl << "PlotTest::";
	// std::cout << "signal_spectr() " << std::endl;
	// PlotTest::signal_spectr();

	std::cout << std::endl << "PlotTest::";
	std::cout << "plot_energy(z=5,10,20) " << std::endl;
	PlotTest::plot_energy(5);
	PlotTest::plot_energy(10);
	PlotTest::plot_energy(20);
}

void PlotTest::signal_spectr ()
{
	double tau = 0.5;
	auto f = [tau] (double vt) { return Function::gauss(vt,tau); };

	for (double w0 = 1e5; w0 < 1e7; w0 += 1) {
		auto re = [f,w0] (double w) { return f(w) * cos(w0 * w); };
		auto im = [f,w0] (double w) { return f(w) * sin(w0 * w); };
		double spec = 0;
		double norm = Math::inv_sqrt(2*M_PI);
		spec += SimpsonRunge(1e2,1,1e4).value(0,1e3,re);
		spec -= SimpsonRunge(1e2,1,1e4).value(0,1e3,im);
		if (spec) std::cout << w0 << " -> " << spec/norm << std::endl;
	}
}

void PlotTest::plot_energy (double z)
{
	double R = 1, A0 = 1, tau = 0.5;
	double eps_r = 1, mu_r = 1;

	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(R, A0);
	MissileField* linear = new MissileField(source, medium);
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source, tau);
	free_shape->set_time_depth([tau] (double vt) {return Function::gauss(vt,tau);});
	LinearDuramel* duhamel = new LinearDuramel(free_shape, medium, linear);

	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	double max  = ((Electrodynamics*) duhamel)->energy_e(0.0,0.0,z);
	std::cout << std::endl;
	std::cout << "E(0,z=" << z << "): " << max << std::endl;

	for (double x = -3; x <= 3;  x += 0.1) {
		for (double y = -3; y <= 3;  y += 0.1) {
			double rho = std::sqrt(x*x + y*y);
			double phi = std::atan2(y, x);
			double val = ((Electrodynamics*) duhamel)->energy_e(rho,phi,z) / max;
			std::cout << x << ' ' << y << " -> " << val << std::endl;
			line = {x, y, val};
			plot_data.push_back(line);
			line.clear();
		}
	}

	GnuPlot* plot = new GnuPlot( std::to_string((int)z) + ".gnp" );
	plot->set_gnuplot_bin( PlotTest::global_conf->path_gnuplot_binary() );
	plot->set_ox_label("x, m");
	plot->set_oy_label("y, m");
	plot->grid_on();
	plot->cage_on();
	plot->plot_colormap(plot_data);
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

void PlotTest::set_options ()
{
	PlotTest::global_conf->path_gnuplot_binary("gnuplot/bin/gnuplot");
	/* TODO: BUG - not used config param
	PlotTest::global_conf->plot_color_map(Colormap::parula); */
}

void PlotTest::plot_source ()
{
	std::vector<std::vector<double>> plot_data;
	std::vector<double> line;

	// auto f = [] (double x) { return Function::rect(x, 1); };
	// auto f = [] (double x) { return Function::sin(x, 1); };
	// auto f = [] (double x) { return Function::sinc(x, 1); };
	// auto f = [] (double x) { return Function::gauss(x, 1); };
	auto f = [] (double x) { return Function::sigmoid(x, 1); };
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
	plot->set_oy_label("f(t)");
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
	plot->set_colormap(Colormap::parula);
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
