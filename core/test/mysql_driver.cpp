//
//  mysql_driver.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 07.03.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "mysql_connect.hpp"
using namespace std;

#define HOSTNAME "localhost"
#define USERNAME "maxwell"
#define PASSWORD "maxwell"
#define DATABASE "maxwell"

#define COMMENT "unit test problem 42"

static size_t INSERTED_POINTS;
static size_t TEST_PROBLEM_ID;

bool create_problem ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    TEST_PROBLEM_ID = client.insert_problem(COMMENT);

    double rho = 0, phi = 0, z = 2;
    for (int ct = 0; ct < 50; ct++) {
        client.select_probe(ct,rho,phi,z);
        client.update_probe_result(ct);
        INSERTED_POINTS++;
    }

    return (client.select_all_probes().size() == INSERTED_POINTS) ? true : false;
}

bool find_problem ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    auto saved = client.get_saved_problems();

    for (auto problem : saved) {
        if (!problem.second.compare(COMMENT)) {
            if (problem.first == TEST_PROBLEM_ID)
                return true;
        }
    }

    return false;
}

bool check_data ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    client.select_problem(TEST_PROBLEM_ID);

    size_t point_number = 0;
    double rho = 0, phi = 0, z = 2;
    for (int ct = 0; ct < 50; ct++) {
        client.select_probe(ct,rho,phi,z);
        if (client.get_probe_result() == ct) {
            point_number++;
        }
    }

    if (client.select_all_probes().size() != INSERTED_POINTS) return false;
    if (point_number != INSERTED_POINTS) return false;
    return true;
}

bool delete_problem ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    client.select_problem(TEST_PROBLEM_ID);

    double rho = 0, phi = 0, z = 2;
    for (int ct = 0; ct < 50; ct++) {
        client.select_probe(ct,rho,phi,z);
        client.delete_probe();
    }

    bool match = (client.select_all_probes().size() == 0) ? true : false;

    client.delete_problem();

    return match;
}

size_t insert_evolution_data ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    TEST_PROBLEM_ID = client.insert_problem(COMMENT);
    client.select_probe(0,0,0,0);

    size_t inserted = 0;
    for (float nu = 0; nu < 50; nu += 0.5) {
        client.select_coeff(-1,nu);
        client.update_coeff_result(nu);
        inserted += 1;
    }

    return inserted;
}

bool ckeck_evolution_data (size_t inserted)
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    client.select_problem(TEST_PROBLEM_ID);
    client.select_probe(0,0,0,0);

    size_t matched = 0;
    for (float nu = 0; nu < 50; nu += 0.5) {
        client.select_coeff(-1,nu);
        if (client.get_coeff_result() == nu) matched++;
        client.delete_coeff();
    }

    bool empty_coeff = (client.select_all_coeffs().first.size() == 0) ? true : false;

    // client.delete_probe(); // Not nessesery
    client.delete_problem();

    return empty_coeff && (matched == inserted);
}

int main ()
{
    if (!create_problem()) return 1;
    if (!find_problem()) return 1;
    if (!check_data()) return 1;
    if (!delete_problem()) return 1;

    auto inserted = insert_evolution_data();
    if (ckeck_evolution_data(inserted)) return 1;

    return 0;
}
