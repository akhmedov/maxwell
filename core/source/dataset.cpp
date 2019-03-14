//
//  dataset.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"

Dataset::Dataset (std::size_t rdx, double dc, double noise_pw)
: radix(rdx), noise_power(noise_pw), duty_cycle(dc), noise(0,std::sqrt(noise_power))
{ }

Dataset::~Dataset ()
{}

double Dataset::count_snr_db (std::vector<double> field, double noise_power)
{
	double power = 0;
	for (auto i : field)
		power += i*i;
	power /= field.size();
	return 10 * log10(power/noise_power);
}

void Dataset::append (const std::vector<std::vector<double>>& arg, const std::vector<double>& field, std::size_t id)
{
	Dataset::Metadata info;
	info.space = std::vector<double>(arg[0].begin()+1, arg[0].end());
	info.signal = id;
	info.snr = Dataset::count_snr_db(field, noise_power);
	script.push_back(info);

	double absolut_from = series.rbegin()->first + arg[0][1] - arg[0][0];
	double relative_from = arg[0][0];

	for (std::size_t i = 0; i < field.size(); i++) {
		SeriesItem item;
		item.magnitude = field[i] + noise.value(0,0,0,0,0);
		item.info.push_back(script.back());
		double absolute_time = absolut_from + arg[i][0] - relative_from;
		series[absolute_time] = item;
	}

	size_t period = (size_t) duty_cycle * field.size();
	this->append_duty_cycle(arg[0][1] - arg[0][0], period - field.size());
}

void Dataset::append_duty_cycle (double dlt_time, std::size_t points)
{
	struct Metadata zero;
	double absolut_from = series.rbegin()->first + dlt_time;

	for (std::size_t i = 0; i < points; i++) {
		SeriesItem item;
		item.magnitude = noise.value(0,0,0,0,0);
		item.info.push_back(zero);
		double absolute_time = absolut_from + dlt_time * i;
		series[absolute_time] = item;
	}
}

void Dataset::superpos (const std::vector<std::vector<double>>& arg, const std::vector<double>& field, std::size_t id, double crossection)
{
	throw std::logic_error("Not implemented!");
}

std::size_t Dataset::get_radix () const
{
	return radix;
}

std::vector<std::size_t> Dataset::get_signal_at (double time) const
{
	std::vector<std::size_t> res;

	for (auto i : series.lower_bound(time)->second.info)
		res.push_back(i.signal);

	return res;
}

double Dataset::get_amplitude_at (double time) const
{
	return series.lower_bound(time)->second.magnitude;
}

nlohmann::json Dataset::to_json (const Metadata& meta)
{
	nlohmann::json js {
		{"unit", meta.signal},
		{"SNR", meta.snr},
		{"space", meta.space}
	};

	return js;
}

nlohmann::json Dataset::to_json (const std::map<double, SeriesItem>& series)
{
	std::vector<nlohmann::json> js;

	for (auto item : series) {

		std::vector<nlohmann::json> meta;
		for (auto info : item.second.info)
			meta.push_back(Dataset::to_json(info));

		nlohmann::json tmp {
			{"time", item.first},
			{"emp", item.second.magnitude},
			{"sparks", item.second.info.size()},
			{"meta", meta}
		};

		js.push_back(tmp);
	}

	return js;
}

nlohmann::json Dataset::get_dataset (const std::string& name) const
{
	nlohmann::json dataset {
		{"name", name},
		{"radix", radix},
		{"duty_cycle", duty_cycle},
		{"noise_pw", noise_power},
		{"sparks", script.size()},
		{"points", series.size()},
		{"series", Dataset::to_json(series)}
	};

	return dataset;
}

void Dataset::serialize (const std::string& filename, const nlohmann::json& dataset, bool binary)
{
	if (binary) {
		// std::vector<std::uint8_t> v_ubjson = json::to_ubjson(j);
		// json j_from_ubjson = json::from_ubjson(v_ubjson);
		throw std::logic_error("binary mode is not implemented in Dataset::serialize");
	}

	std::ofstream file;
	file.open(filename);
	file << std::setw(4) << dataset << std::endl;
	file.close();
}
