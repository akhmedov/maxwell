//
//  mysql_driver.cpp
//  test.core.maxwell
//
//  Created by Rolan Akhmedov on 07.03.19.
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "maxwell.hpp"
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
        client.select_point(ct,rho,phi,z);
        client.update_result(ct);
        INSERTED_POINTS++;
    }

    return (client.select_entity().size() == INSERTED_POINTS) ? true : false;
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
        client.select_point(ct,rho,phi,z);
        if (client.get_result() == ct) {
            point_number++;
        }
    }

    if (client.select_entity().size() != INSERTED_POINTS) return false;
    if (point_number != INSERTED_POINTS) return false;
    return true;
}

bool delete_problem ()
{
    MySQL client(HOSTNAME, USERNAME, PASSWORD, DATABASE);
    client.select_problem(TEST_PROBLEM_ID);

    double rho = 0, phi = 0, z = 2;
    for (int ct = 0; ct < 50; ct++) {
        client.select_point(ct,rho,phi,z);
        client.delete_point();
    }

    bool match = (client.select_entity().size() == 0) ? true : false;

    client.delete_problem();

    return match;
}

int main ()
{
    if (!create_problem()) return 1;
    if (!find_problem()) return 1;
    if (!check_data()) return 1;
    if (!delete_problem()) return 1;

    return 0;
}
