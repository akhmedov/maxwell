//
//  abstract_field.hpp
//  Maxwell
//
//  Created by Oleh Zahrychanskyi on 19.05.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#pragma once

#include <algorithm>
#include <exception>
#include <atomic>
#include <thread>
#include <vector>
#include <cmath>
#include <iostream>

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
// struct DatabaseCalculationManager {

//     DatabaseCalculationManager (std::size_t threads_);
//     ~DatabaseCalculationManager ();

//     void wait ();
//     std::size_t progress () const;

//     template <template <typename ...> class ArgT, typename T, typename FuncT>
//     void start (const std::vector<ArgT<T>>& input, std::vector<T>& output, FuncT func);

//     template <class ArgT, typename ResT, typename ClassT, typename FuncT>
//     void start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property);

// private:

//     template<template<typename ...> class ArgT, typename T, typename FuncT>
//     void startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func);

//     template <class ArgT, typename ResT, typename ClassT, typename FuncT>
//     void startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, std::size_t from, std::size_t to, ClassT* field, FuncT property);

//     std::size_t threads;
//     std::vector<std::thread> thread_pool;
//     std::atomic_size_t current{ 1 };
// };

// ### CalculationManager implementation ##################################################################

inline CalculationManager::CalculationManager (std::size_t threads_)
: threads(threads_) 
{
    thread_pool.reserve(threads_);
}

inline CalculationManager::~CalculationManager ()
{
    wait();
}

inline std::size_t CalculationManager::progress () const 
{
    return current;
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
    if (input.size() < threads) threads = input.size();

    current = 0;
    int from{};
    auto len = input.size() / threads;
    if (input.size() % threads) len++;
    for (auto i = 0u; i < threads; i++) {
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
        current++;
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
        current++;
    }
}

// ### DatabaseCalculationManager implementation ########################################################

// inline DatabaseCalculationManager::DatabaseCalculationManager (std::size_t threads_)
// : threads(threads_) 
// {
//     thread_pool.reserve(threads_);
// }

// inline DatabaseCalculationManager::~DatabaseCalculationManager ()
// {
//     wait();
// }

// inline std::size_t DatabaseCalculationManager::progress () const 
// {
//     return current;
// }

// inline void DatabaseCalculationManager::wait () 
// {
//     for (auto &thread : thread_pool)
//         if (thread.joinable())
//             thread.join();
// }

// template <class ArgT, typename ResT, typename ClassT, typename FuncT>
// void DatabaseCalculationManager::start (const std::vector<ArgT>& input, std::vector<ResT>& output, ClassT* field, FuncT property)
// {
//     if (!field) throw std::logic_error("filed can not be NULL or nullptr");
//     if (input.size() != output.size()) throw std::logic_error("Size of argument and result does not match");
//     if (input.size() < threads) threads = input.size();

//     current = 0;
//     int from{};
//     auto len = input.size() / threads;
//     if (input.size() % threads) len++;
//     for (auto i = 0u; i < threads; i++) {
//         auto to = std::min(from + len, input.size());
//         thread_pool.emplace_back(
//             &DatabaseCalculationManager::startImpl<ArgT, ResT, ClassT, FuncT>, this, 
//             std::ref(input), std::ref(output), from, to, field, property);
//         from = to;
//     }
// }

// template <class ArgT, typename ResT, typename ClassT, typename FuncT>
// void DatabaseCalculationManager::startImpl (const std::vector<ArgT> &input, std::vector<ResT> &output, std::size_t from, std::size_t to, ClassT* field, FuncT property)
// {
//     auto fnc = static_cast<std::function<ResT(ClassT*,ArgT)>>(property);
//     for (auto i = from; i < to; i++) {
//         output[i] = fnc(field,input[i]);
//         current++;
//     }
// }

// template<template<typename ...> class ArgT, typename T, typename FuncT>
// void DatabaseCalculationManager::start (const std::vector<ArgT<T>> &input, std::vector<T> &output, FuncT func)
// {
//     if (input.size() != output.size()) throw std::logic_error("Size of argument and result does not match");
//     if (input.size() < threads) threads = input.size();

//     current = 0;
//     int from{};
//     auto len = input.size() / threads;
//     if (input.size() % threads) len++;
//     for (auto i = 0u; i < threads; i++) {
//         auto to = std::min(from + len, input.size());
//         thread_pool.emplace_back(
//             &DatabaseCalculationManager::startImpl<ArgT, T, FuncT>, this, 
//             std::ref(input), std::ref(output), from, to, func);
//         from = to;
//     }
// }

// template<template<typename ...> class ArgT, typename T, typename FuncT>
// void DatabaseCalculationManager::startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func)
// {
//     for (auto i = from; i < to; i++) {
//         output[i] = func(input[i]);
//         current++;
//     }
// }
