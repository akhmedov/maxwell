//
//  dataset.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 26.06.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef dataset_hpp
#define dataset_hpp

#define RADIUS 1
#define AMPLITUDE 1

#include "noise.hpp"
#include "manager.hpp"
#include "function.hpp"
#include "nlohmann_json.hpp"
#include "abstract_field.hpp"
#include "linear_duhamel.hpp"

#include <algorithm>
#include <random>
#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include <string>

using json = nlohmann::json;

typedef std::pair<double,double> min_max;
typedef std::function<double(AbstractField*,double,double,double,double)> Component;

namespace serial {

	struct dataset {

		dataset (const std::vector<AbstractField*>& emp_shape, 
				 double effective_duration,
				 AdditiveWhiteGaussian* noise = NULL, 
				 double duty_cycle = 0.5, 
				 double vt_step = 0.01);

		void set_comonent (const Component& comp);
		void set_char (double rho, double phi, double z, std::size_t id);
		void set_char (double rho, double phi, double z, std::size_t id1, std::size_t id2, double cross_section);

		double get_real_vt (std::size_t id_vt);
		double get_rho (std::size_t id_vt);
		double get_phi (std::size_t id_vt);
		double get_z (std::size_t id_vt);

		double get_amplitude (std::size_t id_vt);
		std::size_t get_sparks ();
		std::size_t get_radix  ();
		std::vector<std::size_t> get_signal_at (std::size_t id_vt);
		std::vector<double> get_series ();

	protected:

		void evaluate ();
		static AbstractField* updisk_field (const std::function<double(double)>& f);
		static double updisk_emp_duraton (double tau, double rho, double z);

	private:

		std::size_t sparks_number;
		double time_step;
		double duty_cycle;
		double effective_duration;
		Component component;
		AdditiveWhiteGaussian* background;
		std::function<std::vector<double>(double)> route;
		std::vector<AbstractField*> characters;
		std::vector<std::vector<double>> series;
		// series[i] = {vt, field, real_vt, rho, phi, z, id1}
	};

	// void randomized_sequental (std::size_t pulses, std::size_t radix, double sigma, 
	// 							const std::string& file_name,
	// 							min_max rho, min_max phi, min_max z);
	// void same_snr (std::size_t pulses, std::size_t radix, double snr, 
	// 							const std::string& file_name,
	// 							double rho, double phi, double z);
	json json_from (dataset series, min_max rho, min_max phi, min_max z);
	void serialize (const std::string& filename, const json& dataset, bool binary = false);
};

#endif /* dataset_hpp */
