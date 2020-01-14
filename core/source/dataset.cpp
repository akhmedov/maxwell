//
//  dataset.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"

void Dataset::to_json (nlohmann::json& js, const Annotation& data)
{
	js = {
		{"signal_id", data.signal_id},
		{"time_step", data.time_step},
		{"SNR", data.snr_ratio},
		{"observ_point", data.space_point},
		{"time_from", data.time_idx_from},
		{"time_to", data.time_idx_to}
	};
}

void Dataset::to_json (nlohmann::json& js, const SeriesItem& item)
{
	js = {
		{"field", item.field},
		{"annotation", item.annotation}
	};
}

Dataset::Dataset::Dataset (std::size_t rdx, std::size_t item_lth, double noise_pw, const std::vector<std::string>& cl_label)
: radix(rdx), item_length(item_lth), noise_power(noise_pw), noise(0,std::sqrt(noise_pw)), class_label(cl_label)
{
	this->generator = std::mt19937(this->random_device());
}

void Dataset::Dataset::append (const std::vector<double>& field, const Annotation& data)
{
	if (field.size() >= this->item_length)
		std::invalid_argument("Can not append item to dataset: dataset item ");

	std::vector<double> new_field(item_length, 0);
	for (auto& i : new_field) i += noise.value();

	auto rand = std::uniform_int_distribution<std::size_t>(0,this->item_length - field.size() - 1);
	std::size_t start = rand(this->generator);
	
	std::transform(new_field.begin() + start,  new_field.begin() + start + field.size(),  field.begin(), new_field.begin() + start,  std::plus<double>());
	
	this->dataset.emplace_back(new_field, data);
	this->dataset.back().annotation.time_idx_from = start;
	this->dataset.back().annotation.time_idx_to = start + field.size();
	this->dataset.back().annotation.snr_ratio = Dataset::count_snr_db(field, this->noise_power);
}

double Dataset::Dataset::count_snr_db (std::vector<double> field, double noise_power)
{
	double signal_pawer = 0;
	for (auto i : field) signal_pawer += i*i;
	signal_pawer /= field.size();
	return 10 * log10(signal_pawer/noise_power);
}

void Dataset::Dataset::serialize_to_json (std::string filename, bool binary)
{
	nlohmann::json js = {
		{"name", "no_data"},
		{"radix", this->radix},
		{"item_length", this->item_length},
		{"items_in_dataset", this->dataset.size()},
		{"noise_power", this->noise_power},
		{"class_label", this->class_label},
		{"dataset", this->dataset}
	};

	if (binary) {
		// std::vector<std::uint8_t> v_ubjson = json::to_ubjson(j);
		// json j_from_ubjson = json::from_ubjson(v_ubjson);
		throw std::logic_error("binary mode is not implemented in Dataset::serialize");
	}

	std::ofstream file;
	file.open(filename);
	file << std::setw(4) << js << std::endl;
	file.close();	
}

Dataset::Dataset::~Dataset() {}

/* nlohmann::json Dataset::read_file (const std::string& filename)
{
	std::ifstream file(filename);
	return nlohmann::json::parse(file);
} */
