//
//  script_manager.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 16.10.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "script_manager.hpp"

std::vector<std::vector<double>> ScriptManager::datagrid_from (const std::vector<std::pair<double,double>> &args, const std::vector<double> &fnc, double &xmin, double &xstep, double &xmax, double &ymin, double &ystep, double &ymax)
{
	std::map<double, std::map<double,double>> maped; // y -> {x -> F}

	if (fnc.size() != args.size()) 
		throw std::invalid_argument("Size of arg and fnc vectors must matchs.");

	for (size_t i = 0; i < fnc.size(); i++) {
		double x = args[i].first;
		double y = args[i].second;
		maped[y][x] = fnc[i];
	}

	xmin = maped.begin()->second.begin()->first;
	xstep = std::next(maped.rbegin()->second.begin())->first - maped.rbegin()->second.begin()->first;
	xmax = maped.rbegin()->second.rbegin()->first;
	ymin = maped.begin()->first;
	ystep = std::next(maped.begin())->first - maped.begin()->first;
	ymax = maped.rbegin()->first;

	std::vector<std::vector<double>> grid;

	while (!maped.empty()) {

		auto sallestY = maped.begin();
		std::map<double,double> map_sameY = sallestY->second; // get submap with same Y
		maped.erase(sallestY->first); // erase submap from master map

		std::vector<double> vec_sameY; // transform map to vector
		for (auto it = map_sameY.begin(); it != map_sameY.end(); ++it)
			vec_sameY.push_back( it->second );

		grid.push_back(vec_sameY); // add line to result matrix
	}

	size_t items = 0;
	for (auto&& line : grid) items += line.size();

	return grid;
}

void ScriptManager::write_script (const std::string &text) const
{
	std::ofstream script;
	script.open( this->fname );
	script << text;
	script.close();

	if (this->log) {
		std::string msg = "Write comands to $FNAME script... Done.";
		msg = std::regex_replace(msg, std::regex("\\$FNAME"), this->fname);
		this->log->info(msg);
	}
}
