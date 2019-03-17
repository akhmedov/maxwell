//
//  dataset_series.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 15.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"
#include "nlohmann_json.hpp"

#include <string>
#include <numeric>
#include <iostream>

using namespace std;

static const int    RADIX       = 2;
static const int    SPARKS_NUM  = 10;
static const int    POINTS_NUM  = 100; // TODO:
static const double DUTY_SYCLE  = 0.5;
static       double DATASUM     = 0;
static const double NOISE_POWER = 1.0;
static const string DATASETNAME = "dump";
static const string TESTFILE    = "test_dataset.json";

void create_dataset ()
{
    Dataset ds(RADIX, DUTY_SYCLE, NOISE_POWER);

    for (int i = 0; i < SPARKS_NUM; i++) {
        vector<double> time {0, 1, 2, 3, 4, 5, 6};
        vector<double> func {0, 1, 2, 3, 4, 5, 0}; 
        DATASUM += accumulate(func.begin(), func.end(), 0);
        ds.append(1, {M_PI,M_PI,M_PI}, time, func);
    }

    nlohmann::json js = ds.get_dataset(DATASETNAME);
    Dataset::serialize(TESTFILE, js, false);
}

bool check_values (const nlohmann::json& dataset)
{
    double sum = 0;
    double last_time = 0;
    for (auto item : dataset["series"]) {
        if (last_time > item["time"].get<double>()) return 1;
        last_time = item["emp"].get<double>();
        sum += item["emp"].get<double>();
    }
     
    return abs(DATASUM - sum) < 50;
}

bool check_script ()
{
    // TODO: check interal script and parsed script
    return false;
}

bool check_ ()
{
    // TODO: check interal script and parsed script
    return false;
}

int main ()
{
    create_dataset();
    nlohmann::json dataset = Dataset::read_file(TESTFILE);
    if (!check_values(dataset)) return 1;
    return 0;
} 
