//
//  dataset.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"

void DatasetHelper::to_json (nlohmann::json& js, const DatasetHelper::Metadata& meta)
{
	js = {
		{"unit", meta.signal},
		{"SNR", meta.snr},
		{"space", meta.space}
	};
}

void DatasetHelper::to_json (nlohmann::json& js, const DatasetHelper::SeriesItem& item)
{
	js = {
		{"time", 0},
		{"magnitute", item.magnitude},
		{"info", item.info}
	};
}

Dataset::Dataset (std::size_t rdx, double dc, double noise_pw)
: radix(rdx), duty_cycle(dc), noise_power(noise_pw), noise(0,std::sqrt(noise_power))
{
	script.emplace_back();

	for (double time = 0; time < 2; time += 0.01)
		series[time] = {noise.value(), {script.front()}};
}

double Dataset::count_snr_db (std::vector<double> field, double noise_power)
{
	double signal_pawer = 0;
	for (auto i : field) signal_pawer += i*i;
	signal_pawer /= field.size();
	return 10 * log10(signal_pawer/noise_power);
}

void Dataset::append (std::size_t id, const std::vector<double>& space, const std::vector<double>& time, const std::vector<double>& field)
{
	if (id > radix - 1)
		throw std::logic_error("Not legal id");
	if (time.size() != field.size())
		throw std::logic_error("Size of time and field are not matched");
	if (time.size() < 2)
		throw std::logic_error("Data size to low");

	double snr = Dataset::count_snr_db(field, noise_power);
	script.emplace_back(snr, id, space);

	double time_dlt = std::abs(time[1] - time[0]);
	double absolut = series.rbegin()->first + time_dlt;
	double relative = time[0];

	for (std::size_t i = 0; i < field.size(); i++) {
		DatasetHelper::SeriesItem item;
		item.magnitude = field[i] + noise.value();
		item.info.push_back(script.back());
		series[absolut + time[i] - relative] = item;
	}

	std::size_t wait_duration = field.size() / duty_cycle - field.size();
	this->append_duty_cycle(time_dlt, wait_duration);
}

void Dataset::append_duty_cycle (double dlt_time, std::size_t points)
{
	double absolut_from = series.rbegin()->first + dlt_time;

	for (std::size_t i = 0; i < points; i++) {
		DatasetHelper::SeriesItem item;
		item.magnitude = noise.value();
		item.info.emplace_back();
		series.emplace(absolut_from + dlt_time * i, std::move(item));
	}
}

void Dataset::superpos (std::size_t /*id*/, const std::vector<double>& /*space*/, const std::vector<double>& /*time*/, const std::vector<double>& /*field*/, double /*crossection*/)
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

nlohmann::json Dataset::to_json (const std::map<double, DatasetHelper::SeriesItem>& series)
{
	std::vector<nlohmann::json> js;

	for (auto item : series) {

		std::vector<nlohmann::json> meta;
		for (auto info : item.second.info)
			meta.push_back(info);

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
	nlohmann::json dataset = {
		{"name", name},
		{"radix", radix},
		{"duty_cycle", duty_cycle},
		{"noise_pw", noise_power},
		{"sparks", script.size()-1},
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

nlohmann::json Dataset::read_file (const std::string& filename)
{
	std::ifstream file(filename);
	return nlohmann::json::parse(file);
}

Dataset* Dataset::instance_from (nlohmann::json input)
{
	UNUSED(input);
	throw std::logic_error("Not implemented!");

	Dataset* output = new Dataset(0,0,0);
	return output;
}
