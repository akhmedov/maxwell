#pragma once

#include <algorithm>
#include <exception>
#include <atomic>
#include <thread>
#include <vector>
#include <cmath>

// NOTE: template functions have `typename ...` paramiter to be comartable 
// with STL containers like vector and custom classes like Data

// Oneshot cluster to calculate data in parallel.
struct CalculationManager {

    CalculationManager (std::size_t threads_);
    ~CalculationManager ();

    void wait ();
    std::size_t progress () const;

    template<template<typename ...> class ArgT, typename T, typename FuncT>
    void start (const std::vector<ArgT<T>> &input, std::vector<T> &output, FuncT func);

private:

    template<template<typename ...> class ArgT, typename T, typename FuncT>
    void startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func);

    std::size_t threads;
    std::vector<std::thread> thread_pool;
    std::atomic_size_t current{ 1 };
};

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
