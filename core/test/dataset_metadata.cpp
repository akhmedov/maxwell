//
//  dataset_metadata.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 15.03.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "dataset.hpp"
#include "nlohmann_json.hpp"

#include <string>
#include <numeric>

using namespace std;

static const int    RADIX       = 2;
static const int    SPARKS_NUM  = 10;
// static const int    POINTS_NUM  = 100; // TODO: add testing
static const double DUTY_SYCLE  = 0.5;
static       double DATASUM     = 0;
static const double NOISE_POWER = 1.0;
static const string DATASETNAME = "dump";
static const string TESTFILE    = "test_dataset.json";

bool floating_compare (double f1, double f2)
{
    return (abs(f1-f2) / f2) < 0.01;
}

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

bool check_metadata (const nlohmann::json& dataset)
{
    bool nm_comp = !dataset["name"].get<string>().compare(DATASETNAME);
    bool dc_comp = floating_compare(dataset["duty_cycle"], DUTY_SYCLE);
    bool np_comp = floating_compare(dataset["noise_pw"], NOISE_POWER);
    bool nw_comp = NOISE_POWER == dataset["noise_pw"].get<int>();
    bool rd_comp = RADIX == dataset["radix"].get<int>();
    bool sp_comp = SPARKS_NUM == dataset["sparks"].get<int>();
    // bool pn_comp = POINTS_NUM == dataset["points"].get<int>();
    return nm_comp && dc_comp && np_comp && nw_comp && rd_comp && sp_comp;
}


int main ()
{
    create_dataset();
    nlohmann::json dataset = Dataset::read_file(TESTFILE);
    if (!check_metadata(dataset)) return 1;
    return 0;
} 
