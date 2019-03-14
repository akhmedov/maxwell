//
//  dataset.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef dataset_hpp
#define dataset_hpp

#include "maxwell.hpp"
#include "nlohmann_json.hpp"

#include <random>
#include <algorithm>
#include <functional>

#include <vector>
#include <list>

#include <fstream> // ofstream
#include <string> // to_string
#include <memory> // unique_ptr

// using json = nlohmann::json;
// using func = std::function<double(std::vector<double>)>;
// using comp = std::function<double(AbstractField*,double,double,double,double)>;

struct Dataset {

	Dataset (std::size_t radix, double duty_cycle, double noise_power);
	~Dataset ();

	void append   (const std::vector<std::vector<double>>& arg, const std::vector<double>& field, std::size_t id);
	void superpos (const std::vector<std::vector<double>>& arg, const std::vector<double>& field, std::size_t id, double crossection);

	// dataset global metadata getters
	std::size_t get_radix () const;
	nlohmann::json get_dataset (const std::string& name) const;

	// point getters
	std::vector<std::size_t> get_signal_at (double time) const;
	double get_amplitude_at (double time) const;

	static void serialize (const std::string& filename, const nlohmann::json& dataset, bool binary);

private:

	struct Metadata {
		double snr {0};
		std::size_t signal {0};
		std::vector<double> space {0,0,0};
	};

	struct SeriesItem {
		double magnitude {};
		std::vector<Metadata> info {}; // TODO: not a copy but ref to script[i]
	};

	void append_duty_cycle (double dlt_time, std::size_t points);
	static nlohmann::json to_json (const Metadata& meta);
	static nlohmann::json to_json (const std::map<double, SeriesItem>& series);
	static double count_snr_db (std::vector<double> filed, double noise_power);

	std::size_t radix;
	double duty_cycle;
	double noise_power;
	AdditiveWhiteGaussian noise;
	std::vector<Metadata> script;
	std::map<double, SeriesItem> series;
};

#endif /* dataset_hpp */
