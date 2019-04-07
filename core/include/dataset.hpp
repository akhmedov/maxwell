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

#include <list>
#include <vector>
#include <fstream> // ofstream

namespace DatasetHelper {

	struct Metadata {
		Metadata () : Metadata (0,0,{NAN,NAN,NAN}) {}
		Metadata (double ratio, std::size_t id, std::vector<double> coord) : snr(ratio), signal(id), space(coord) {}
		double snr {0};
		std::size_t signal {0};
		std::vector<double> space {0,0,0};
	};

	void to_json (nlohmann::json& js, const Metadata& meta);

	struct SeriesItem {
		double magnitude {};
		std::vector<Metadata> info {}; // TODO: not a copy but pointer to inctance on heap
	};

	void to_json (nlohmann::json& js, const SeriesItem& item);
};

struct Dataset {

	Dataset (std::size_t radix, double duty_cycle, double noise_power);

	void append   (std::size_t id, const std::vector<double>& space, const std::vector<double>& time, const std::vector<double>& field);
	void superpos (std::size_t id, const std::vector<double>& space, const std::vector<double>& time, const std::vector<double>& field, double crossection);

	// dataset global metadata getters
	std::size_t get_radix () const;
	nlohmann::json get_dataset (const std::string& name) const;

	// point getters
	std::vector<std::size_t> get_signal_at (double time) const;
	double get_amplitude_at (double time) const;

	static void serialize (const std::string& filename, const nlohmann::json& dataset, bool binary);
	static nlohmann::json read_file (const std::string& filename);
	static Dataset* instance_from (nlohmann::json);

private:

	void append_duty_cycle (double dlt_time, std::size_t points);
	static nlohmann::json to_json (const std::map<double, DatasetHelper::SeriesItem>& series);
	static double count_snr_db (std::vector<double> filed, double noise_power);

	std::size_t radix;
	double duty_cycle;
	double noise_power;
	WhiteGaussian noise;
	std::list<DatasetHelper::Metadata> script; // signals cannt come at the same time
 	std::map<double, DatasetHelper::SeriesItem> series;
};

#endif /* dataset_hpp */
