//
// Created by rui on 11/16/18.
//

#ifndef STAT906A3_GLOBAL_H
#define STAT906A3_GLOBAL_H

#include "matplotlibcpp.h"
#include <random>
namespace plt = matplotlibcpp;

std::random_device rd;
std::mt19937 rgen(rd());
using std::normal_distribution;
using std::uniform_real_distribution;

#endif //STAT906A3_GLOBAL_H
