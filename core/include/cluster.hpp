#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>
#include <cmath>

// Oneshot cluster to calculate data in parallel.
class Cluster {
public:
    Cluster (std::size_t threads_);
    ~Cluster ();

    void wait ();
    std::size_t progress () const;

    template<template<typename> class ArgT, typename T, typename FuncT>
    void start (const std::vector<ArgT<T>> &input, std::vector<T> &output, FuncT func);

private:

    template<template<typename> class ArgT, typename T, typename FuncT>
    void startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func);

    std::size_t threads;
    std::vector<std::thread> thread_pool;
    std::atomic_size_t current{ 1 };
};

inline Cluster::Cluster (std::size_t threads_)
: threads(threads_) 
{
    thread_pool.reserve(threads_);
}

inline Cluster::~Cluster ()
{
    wait();
}

inline std::size_t Cluster::progress () const 
{
    return current;
}

inline void Cluster::wait () 
{
    for (auto &thread : thread_pool)
        if (thread.joinable())
            thread.join();
}

template<template<typename> class ArgT, typename T, typename FuncT>
void Cluster::start (const std::vector<ArgT<T>> &input, std::vector<T> &output, FuncT func)
{
    // TODO: Check for sizes of input and output.
    current = 0;
    int from{};
    auto len = input.size() / threads;
    if (input.size() % threads) len++;
    for (auto i = 0u; i < threads; i++) { // TODO: What if number of threads is less than input data.
        auto to = std::min(from + len, input.size());
        thread_pool.emplace_back(
            &Cluster::startImpl<ArgT, T, FuncT>, this, 
            std::ref(input), std::ref(output), from, to, func);
        from = to;
    }
}

template<template<typename> class ArgT, typename T, typename FuncT>
void Cluster::startImpl (const std::vector<ArgT<T>> &input, std::vector<T> &output, std::size_t from, std::size_t to, FuncT func)
{
    for (auto i = from; i < to; i++) {
        output[i] = func(input[i]);
        current++;
    }
}
