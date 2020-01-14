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
#include <random>
#include <iostream>
#include <fstream> // ofstream
#include <algorithm> // transform plus

namespace Dataset {

	struct Annotation {
		Annotation (std::size_t id, double time_dlt, std::vector<double> coord) 
			: snr_ratio(0), time_step(time_dlt), signal_id(id), space_point(coord), time_idx_from(0), time_idx_to(0) {}
		double snr_ratio {0};
		double time_step {0};
		std::size_t signal_id {0};
		std::vector<double> space_point {0,0,0};
		std::size_t time_idx_from {0};
		std::size_t time_idx_to {0};
	};

	void to_json (nlohmann::json& js, const Annotation& data);

	struct SeriesItem {
		SeriesItem(const std::vector<double>& f, const Annotation& info)
		: field(f), annotation(info) { }
		std::vector<double> field;
		Annotation annotation;
	};

	void to_json (nlohmann::json& js, const SeriesItem& data);

	struct Dataset {

		Dataset (std::size_t radix, std::size_t item_length, double noise_power, const std::vector<std::string>& cl_label);
		void append (const std::vector<double>& field, const Annotation& data);
		void serialize_to_json (std::string file_name, bool binary = false);
		~Dataset();

		static double count_snr_db (std::vector<double> field, double noise_power);

	private:

		std::random_device random_device;
    	std::mt19937 generator;
		
		std::size_t radix;
		std::size_t item_length;
		double noise_power;
		WhiteGaussian noise;
		std::vector<std::string> class_label;
		std::vector<SeriesItem> dataset;
	};

}

#endif /* dataset_hpp */
