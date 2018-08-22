//
//  dataset.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"

serial::dataset::dataset (const std::vector<std::function<double(double)>>& emp_shape, 
						double ed, AdditiveWhiteGaussian* noise, double dc, double step)
: background(noise) {

	this->effective_duration = ed;
	this->time_step = step;
	this->component = &AbstractField::electric_x;
	this->characters = {new ZeroField()};
	this->duty_cycle = dc;
	this->sparks_number = 0;

	for (auto i : emp_shape)
		this->characters.push_back(dataset::updisk_field(i));

	double total_period = dataset::updisk_emp_duraton(this->effective_duration, 1, 1) / this->duty_cycle;
	
	for (double vt = 0; vt <= total_period; vt += step) {
		double noise = (this->background) ? this->background->value(vt, 0, 0, 0) : 0 ;
		this->series.push_back({vt, noise, 0, 0, 0, 0, 0});
	}
}

void serial::dataset::set_comonent (const Component& comp)
{
	this->component = comp;
}

double serial::dataset::updisk_emp_duraton (double tau0, double rho, double z)
{
	double plus = (rho + RADIUS) * (rho + RADIUS) + z*z;
	double minus = (rho - RADIUS) * (rho - RADIUS) + z*z;
	double tau = tau0 + std::sqrt(plus);
	if (rho < RADIUS) return tau - std::abs(z);
	else return tau - std::sqrt(minus);
}

AbstractField* serial::dataset::updisk_field (const std::function<double(double)>& f) 
{
	double eps_r = 1, mu_r = 1;
	Homogeneous* medium = new Homogeneous(mu_r, eps_r);
	UniformPlainDisk* source = new UniformPlainDisk(RADIUS, AMPLITUDE);
	MissileField* on = new MissileField(source, medium);
	auto property = &AbstractField::electric_x;
	FreeTimeCurrent* free_shape = new FreeTimeCurrent(source); 
	free_shape->set_time_depth(f);	
	LinearDuhamel* duhamel = new LinearDuhamel(free_shape, medium, on, NULL);
	return (AbstractField*)duhamel;
}

void serial::dataset::set_char (double rho, double phi, double z, std::size_t id)
{
	if (id >= this->characters.size())
		throw std::logic_error("Signal id is out of characters in signal::set_char()");

	std::uniform_real_distribution<double> uniform(2, 4);
	std::mt19937_64 gen(std::random_device{}());
	double delay = uniform(gen);

	double vt  = this->time_step + this->series.back()[0];
	double vt0 = (rho > RADIUS) ? std::sqrt((rho-RADIUS)*(rho-RADIUS) + z*z) : z;
	double pulse_wides = vt + dataset::updisk_emp_duraton(this->effective_duration, rho, z);
	double total_period = delay + vt + dataset::updisk_emp_duraton(this->effective_duration, rho, z) / this->duty_cycle;

	while (vt < pulse_wides) {
		double noise = (this->background) ? this->background->value(vt, 0, 0, 0) : 0;
		std::vector<double> item = {vt, noise, vt0, rho, phi, z, 0, (double)id};
		this->series.push_back(item);
		vt += this->time_step; vt0 += this->time_step;
	}

	while (vt < total_period) {
		double noise = (background) ? background->value(vt, 0, 0, 0) : 0;
		std::vector<double> item = {vt, noise, 0, 0, 0, 0, 0};
		this->series.push_back(item);
		vt += this->time_step;
	}

	this->sparks_number++;
}

void serial::dataset::set_char (double rho, double phi, double z, std::size_t id1, std::size_t id2, double cross_section)
{
	throw std::logic_error("serial::dataset::set_char is not implemented.");
}

void serial::dataset::evaluate ()
{
	static bool called = false;

	if (!called) {

		std::vector<Manager<0>*> manager;
		for (std::size_t c = 1; c < this->get_radix(); c++) {
			Manager<0>* tmp = new Manager<0>();
			tmp->progress_bar(true);
			manager.push_back(tmp);
		}

		for (std::size_t idx = 0; idx < this->series.size(); idx++) {
			std::vector<std::size_t> unit = this->get_signal_at(idx);
			for (std::size_t c = 1; c < unit.size(); c++) {
				double vt  = this->get_real_vt(idx);
				double rho = this->get_rho(idx);
				double phi = M_PI * this->get_phi(idx) / 180;
				double z   = this->get_z(idx);
				manager[unit[c]-1]->add_argument({vt,rho,phi,z,(double)idx});			
			}
		}

		for (std::size_t c = 1; c < this->get_radix(); c++) {
			auto pair = std::make_pair(this->component,this->characters[c]);
			manager[c-1]->call({pair});
			auto res = manager[c-1]->get_value();
			for (auto i : res) {
				double id    = i[4];
				double field = i[5];
				this->series[id][1] += field;
			}
		}

		called = true;
	}
}

double serial::dataset::get_real_vt (std::size_t vt)
{
	return series[vt][2];
}

double serial::dataset::get_rho (std::size_t vt)
{
	return series[vt][3];
}

double serial::dataset::get_phi (std::size_t vt)
{
	return series[vt][4];
}

double serial::dataset::get_z (std::size_t vt)
{
	return series[vt][5];
}

double serial::dataset::get_amplitude (std::size_t vt)
{
	return series[vt][1];
}

std::size_t serial::dataset::get_sparks ()
{
	return this->sparks_number;
}

std::size_t serial::dataset::get_radix ()
{
	return this->characters.size();
}

std::vector<std::size_t> serial::dataset::get_signal_at (std::size_t vt)
{
	std::vector<double> vd = this->series[vt];
	std::vector<std::size_t> vi = std::vector<std::size_t>(vd.begin(),vd.end());
	vi.erase(vi.begin(),vi.begin()+6);
	return vi;
}

std::vector<double> serial::dataset::get_series ()
{
	this->evaluate();

	std::vector<double> time;
	for (auto i : series)
		time.push_back(i[0]);
	return time;
}

/*

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_int_distribution<> dis(1, 6);
dis(gen) // generation

*/

void serial::randomized_sequental (std::size_t pulses, std::size_t radix, double sigma, const std::string& file_name,
								   min_max rho, min_max phi, min_max z)
{
	double mu = 0;
	double tau = 0.5;
	double duty_cycle = 0.6; 
	double step = 0.01;

	std::vector<std::function<double(double)>> domain;
	domain.push_back([tau] (double vt) { return Function::gauss(vt,tau); });
	// domain.push_back([tau] (double vt) { return Function::sinc(vt,tau); });

	std::mt19937_64 gen1(std::random_device{}());
	std::mt19937_64 gen2(std::random_device{}());
	std::mt19937_64 gen3(std::random_device{}());
	std::mt19937_64 gen4(std::random_device{}());
	std::uniform_real_distribution<double> dist_rho(rho.first, rho.second);
	std::uniform_real_distribution<double> dist_phi(phi.first, phi.second);
	std::uniform_real_distribution<double> dist_z(z.first, z.second);
	std::uniform_int_distribution<>  dist_id(1, domain.size());

	AdditiveWhiteGaussian* noise = new AdditiveWhiteGaussian(mu,sigma);

	serial::dataset series(domain, tau, noise, duty_cycle, step);
	
	for (std::size_t i = 0; i < pulses; i++) {
		std::size_t rnd_id = dist_id(gen1);
		double rnd_rho = dist_rho(gen2);
		double rnd_phi = dist_phi(gen3);
		double rnd_z   = dist_z(gen4);
		series.set_char(rnd_rho, rnd_phi, rnd_z, rnd_id);
	}

	json dataset = serial::json_from(series,rho,phi,z);
	serial::serialize(file_name,dataset);
}

void serial::same_snr (std::size_t pulses, std::size_t radix, double snr, 
	const std::string& file_name, double rho, double phi, double z)
{
	double sigma = 0, Ar = 1, As = 1, Ag = 2;
	double tau = 0.5;
	double duty_cycle = 0.6; 
	double step = 0.01;

	std::vector<std::function<double(double)>> domain;
	domain.push_back([tau,Ag] (double vt) { return Ag * Function::gauss(vt,tau); });
	domain.push_back([tau,As] (double vt) { return As * Function::sinc(vt,tau); });

	std::mt19937_64 gen1(std::random_device{}());
	std::mt19937_64 gen2(std::random_device{}());
	std::mt19937_64 gen3(std::random_device{}());
	std::mt19937_64 gen4(std::random_device{}());
	std::uniform_real_distribution<double> dist_rho(rho, rho);
	std::uniform_real_distribution<double> dist_phi(phi, phi);
	std::uniform_real_distribution<double> dist_z(z, z);
	std::uniform_int_distribution<> dist_id(1, domain.size());

	AdditiveWhiteGaussian* noise = new AdditiveWhiteGaussian(0,sigma);

	serial::dataset series(domain, tau, noise, duty_cycle, step);
	
	for (std::size_t i = 0; i < pulses; i++) {
		std::size_t rnd_id = dist_id(gen1);
		double rnd_rho = dist_rho(gen2);
		double rnd_phi = dist_phi(gen3);
		double rnd_z   = dist_z(gen4);
		series.set_char(rnd_rho, rnd_phi, rnd_z, rnd_id);
	}

	json dataset = serial::json_from(series,std::make_pair(rho,rho),std::make_pair(phi,phi),std::make_pair(z,z));
	serial::serialize(file_name,dataset);
}

json serial::json_from (dataset ds, min_max rho, min_max phi, min_max z)
{
	std::vector<double> time = ds.get_series();
	std::vector<json> series;
	for (std::size_t i = 0; i < time.size(); i++) {
		series.push_back(json{
			{"id",    i},
			{"vt",    time[i]},
			{"field", ds.get_amplitude(i)},
			{"rho",   ds.get_rho(i)},
			{"phi",   ds.get_phi(i)},
			{"z",     ds.get_z(i)},
			{"unit",  ds.get_signal_at(i)}
		});
	}

	json js;
	js["component"] = "unknown";
	js["radix"]     = ds.get_radix();
	js["sparks"]	= ds.get_sparks();
	js["SNR"]       = NAN;
	js["max_rho"]   = rho.second;
	js["min_rho"]   = rho.first;
	js["max_phi"]   = phi.second;
	js["min_phi"]   = phi.first;
	js["max_z"]     = z.second;
	js["min_z"]     = z.first;
	js["series"]    = series;

	return js;
}

void serial::serialize (const std::string& filename, const json& dataset, bool binary)
{
	if (binary) throw std::logic_error("binary mode is not implemented in serial::serialize");
	// std::vector<std::uint8_t> v_ubjson = json::to_ubjson(j);
	// json j_from_ubjson = json::from_ubjson(v_ubjson);
	std::ofstream file;
	file.open(filename);
	file << std::setw(4) << dataset << std::endl;
	file.close();
}
