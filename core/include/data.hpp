#pragma once

#include <vector>

template<typename T> 
struct Space {
    Space() = default;
    Space(T q1_, T q2_, T q3_)
      : q1(q1_), q2(q2_), q3(q3_) {
    }

    T q1{};
    T q2{};
    T q3{};
};

template<typename T>
struct SpaceInterval {
    SpaceInterval() = default;
    // SpaceTime(const)
    SpaceInterval(T q1_, T q2_, T q3_, T from_, T to_)
      : q1(q1_), q2(q2_), q3(q3_), from(from_), to(to_) {
    }

    T q1{};
    T q2{};
    T q3{};
    T from{};
    T to{};
};

using SpaceInterval64 = SpaceInterval<double>;