//
//  abstract_field.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include "mysql_connect.hpp"
#include "space_point.hpp"

#include <typeinfo>
#include <algorithm>
#include <exception>
#include <atomic>
#include <thread>
#include <vector>
#include <cmath>

#include <chrono>
#include <thread>

// NOTE: template functions have `typename ...` paramiter to be comartable 
// with STL containers like vector and custom classes like Data

// #### CalculationManager declaration #####################################################################

// Oneshot cluster to calculate data in parallel.
struct CalculationManager {

    CalculationManager (std::size_t threads_);
    ~CalculationManager ();

    void wait ();
    std::size_t progress () const;

    template <template <typename ...> class ArgT, typename T, typename FuncT>
    void start (const std::vector<ArgT<T>>& input, std::vector<T>& output, FuncT func);

    template <class ArgT, typename ResT, typename ClassT, typename FuncT>
    void start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property);

private:

    template<template<typename ...> class ArgT, typename T, typename FuncT>
    void startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func);

    template <class ArgT, typename ResT, typename ClassT, typename FuncT>
    void startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, std::size_t from, std::size_t to, ClassT* field, FuncT property);

    std::size_t threads;
    std::vector<std::thread> thread_pool;
    std::atomic_size_t current{ 1 };
};

// ### DatabaseCalculationManager declaration #######################################################

// Oneshot cluster to calculate data in parallel with paralel saving to dataset
struct DatabaseCalculationManager {

    DatabaseCalculationManager (std::size_t threads_,
                                std::size_t problem_id,
                                const std::string& server = "localhost",
                                const std::string& database = "maxwell",
                                const std::string& username = "maxwell",
                                const std::string& password = "maxwell");

    DatabaseCalculationManager (std::size_t threads_,
                                const std::string& problem_comment,
                                const std::string& server = "localhost",
                                const std::string& database = "maxwell",
                                const std::string& username = "maxwell",
                                const std::string& password = "maxwell");
    ~DatabaseCalculationManager ();

    void wait ();
    std::size_t progress () const;

    template <template <typename ...> class ArgT, typename T, typename FuncT>
    void start (const std::vector<ArgT<T>>& input, std::vector<T>& output, FuncT func);

    template <class ArgT, typename ResT, typename ClassT, typename FuncT>
    void start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property);

private:

    template <class ArgT, typename ResT, typename ClassT, typename FuncT>
    void startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, const std::vector<std::size_t>& working_idx, ClassT* field, FuncT property, std::size_t thread_id);

    std::size_t threads;
    std::vector<MySQL> connection;
    std::vector<std::thread> thread_pool;
    std::atomic_size_t current{ 1 };
};

// ### CalculationManager implementation ##################################################################

inline CalculationManager::CalculationManager (std::size_t threads_)
: threads(threads_) 
{
    thread_pool.reserve(threads_);
}

inline CalculationManager::~CalculationManager ()
{
    this->wait();
}

inline std::size_t CalculationManager::progress () const 
{
    return this->current;
}

inline void CalculationManager::wait () 
{
    for (auto &thread : thread_pool)
        if (thread.joinable())
            thread.join();
}

template <class ArgT, typename ResT, typename ClassT, typename FuncT>
void CalculationManager::start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property)
{
    if (!field) throw std::logic_error("filed can not be NULL or nullptr");
    if (input.size() != output.size()) throw std::logic_error("Size of argument and result does not match");
    if (input.size() < this->threads) this->threads = input.size();

    current = 0;
    int from{};
    auto len = input.size() / this->threads;
    if (input.size() % this->threads) len++;
    for (auto i = 0u; i < this->threads; i++) {
        auto to = std::min(from + len, input.size());
        thread_pool.emplace_back(
            &CalculationManager::startImpl<ArgT, ResT, ClassT, FuncT>, this, 
            std::ref(input), std::ref(output), from, to, field, property);
        from = to;
    }
}

template <class ArgT, typename ResT, typename ClassT, typename FuncT>
void CalculationManager::startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, std::size_t from, std::size_t to, ClassT* field, FuncT property)
{
    auto fnc = static_cast<std::function<ResT(ClassT*,ArgT)>>(property);
    for (auto i = from; i < to; i++) {
        output[i] = fnc(field,input[i]);
        this->current++;
    }
}

template<template<typename ...> class ArgT, typename T, typename FuncT>
void CalculationManager::start (const std::vector<ArgT<T>> &input, std::vector<T> &output, FuncT func)
{
    if (input.size() != output.size()) throw std::logic_error("Size of argument and result does not match");
    if (input.size() < threads) threads = input.size();

    current = 0;
    int from{};
    auto len = input.size() / threads;
    if (input.size() % threads) len++;
    for (auto i = 0u; i < threads; i++) {
        auto to = std::min(from + len, input.size());
        thread_pool.emplace_back(
            &CalculationManager::startImpl<ArgT, T, FuncT>, this, 
            std::ref(input), std::ref(output), from, to, func);
        from = to;
    }
}

template<template<typename ...> class ArgT, typename T, typename FuncT>
void CalculationManager::startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func)
{
    for (auto i = from; i < to; i++) {
        output[i] = func(input[i]);
        this->current++;
    }
}

// ### DatabaseCalculationManager implementation ########################################################

inline DatabaseCalculationManager::DatabaseCalculationManager (std::size_t threads_,  
                                                                std::size_t problem_id, 
                                                                const std::string& host, 
                                                                const std::string& base, 
                                                                const std::string& user, 
                                                                const std::string& pass)
: threads(threads_)
{
    this->thread_pool.reserve(threads_);

    this->connection.reserve(threads_);
    for (std::size_t i = 0; i < this->threads; i++) {
        this->connection.emplace_back(host, user, pass, base);
        this->connection[i].select_problem(problem_id);
    }
}

inline DatabaseCalculationManager::DatabaseCalculationManager (std::size_t threads_,  
                                                                const std::string& problem_comment, 
                                                                const std::string& host, 
                                                                const std::string& base, 
                                                                const std::string& user, 
                                                                const std::string& pass)
: threads(threads_)
{
    this->thread_pool.reserve(threads_);

    this->connection.reserve(threads_);
    for (std::size_t i = 0; i < this->threads; i++) {
        this->connection.emplace_back(host, user, pass, base);
        this->connection[i].select_problem(problem_comment);
        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
}

inline DatabaseCalculationManager::~DatabaseCalculationManager ()
{
    this->wait();
}

inline std::size_t DatabaseCalculationManager::progress () const 
{
    return this->current;
}

inline void DatabaseCalculationManager::wait () 
{
    for (auto &thread : thread_pool)
        if (thread.joinable())
            thread.join();
}

template <class ArgT, typename ResT, typename ClassT, typename FuncT>
void DatabaseCalculationManager::start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property)
{
    if (!field) throw std::logic_error("filed can not be NULL or nullptr");
    if (input.size() != output.size()) throw std::logic_error("Size of argument and result does not match");
    if (input.size() < this->threads) this->threads = input.size();

    std::vector<std::vector<std::size_t>> working_idx(this->threads, std::vector<std::size_t>());
    std::size_t work = 0;
    for (std::size_t idx = 0; idx < input.size(); idx++) {

        /* if (typeid(input[idx]) == typeid(Point::SpaceTime<Point::Cylindrical>)) {

            this->connection[0].select_probe(input[idx].ct(), input[idx].rho(), input[idx].phi(), input[idx].z());
            double res = this->connection[0].get_probe_result();
            if (!std::isnan(res)) { this->current++; continue; }

        } else if (typeid(input[idx]) == typeid(Point::ModalSpaceTime<Point::Cylindrical>)) {

            this->connection[0].select_probe(input[idx].ct(), input[idx].rho(), input[idx].phi(), input[idx].z());
            this->connection[0].select_coeff(input[idx].m(), input[idx].nu());
            double res = this->connection[0].get_coeff_result();
            if (!std::isnan(res)) { this->current++; continue; }
            
        } else {
            throw std::logic_error("Unknown input item datatype used in DatabaseCalculationManager::start");
        } */

        working_idx[work].push_back(idx);
        work = (work < this->threads-1) ? (work + 1) : 0;
    }

    for (std::size_t i = 0; i < this->threads; i++) {
        if (working_idx[i].empty()) continue;
        thread_pool.emplace_back(
            &DatabaseCalculationManager::startImpl<ArgT, ResT, ClassT, FuncT>, this, 
            std::ref(input), std::ref(output), working_idx[i], field, property, i);
    }
}

template <class ArgT, typename ResT, typename ClassT, typename FuncT>
void DatabaseCalculationManager::startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, const std::vector<std::size_t>& working_idx, ClassT* field, FuncT property, std::size_t thread_id)
{
    auto fnc = static_cast<std::function<ResT(ClassT*,ArgT)>>(property);
    for (auto idx : working_idx) {

        if (typeid(input[idx]) == typeid(Point::SpaceTime<Point::Cylindrical>)) {

            this->connection[thread_id].select_probe(input[idx].ct(), input[idx].rho(), input[idx].phi(), input[idx].z());
            output[idx] = this->connection[thread_id].get_probe_result();
            if (std::isnan(output[idx])) {
                output[idx] = fnc(field,input[idx]);
                this->connection[thread_id].update_probe_result(output[idx]);
            }

        } else if (typeid(input[idx]) == typeid(Point::ModalSpaceTime<Point::Cylindrical>)) {

            this->connection[thread_id].select_probe(input[idx].ct(), input[idx].rho(), input[idx].phi(), input[idx].z());
            this->connection[thread_id].select_coeff(input[idx].m(), input[idx].nu());
            output[idx] = this->connection[thread_id].get_coeff_result();
            if (std::isnan(output[idx])) {
                output[idx] = fnc(field,input[idx]);
                this->connection[thread_id].update_coeff_result(output[idx]);
            }
            
        } else {
            throw std::logic_error("Unknown input item datatype used in DatabaseCalculationManager::start");
        }

        this->current++;
    }
}
