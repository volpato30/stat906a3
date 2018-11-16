//
// Created by rui on 11/16/18.
//

#ifndef STAT906A3_Q3_H
#define STAT906A3_Q3_H
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include "global.h"

class BaseSDE {
public:
    virtual double mu(double x) = 0;
    virtual double sigma(double x) = 0;
    virtual double d_sigma_wrt(double x) = 0;
    virtual ~BaseSDE() = default;;
};

class ShortInterestModel:BaseSDE {
public:
    ShortInterestModel(double kappa, double b, double std, double gamma) {
        _kappa = kappa;
        _b = b;
        _std = std;
        _gamma = gamma;
    }

    double mu(double x) override {
        return _kappa * (_b - x);
    }

    double sigma(double x) override {
        return _std * pow(x, _gamma);
    }

    double d_sigma_wrt(double x) override {
        return _std * _gamma * pow(x, _gamma - 1);
    }

private:
    double _kappa, _b, _std, _gamma;
};

class BaseDiscretize {
public:
    explicit BaseDiscretize(double delta_t, BaseSDE *m) {
        _delta_t = delta_t;
        _normal_gen = normal_distribution<> (0., 1.);
        _model = m;
    }
    virtual double next_x(double x) = 0;

    double get_x_t(double t, double initial_x) {
        auto steps = static_cast<unsigned long>(round(t / _delta_t));
        double temp_x = initial_x;
        for (int i=0; i < steps; i++) {
            temp_x = next_x(temp_x);
        }
        return temp_x;
    }

    double getDelta_t() const {
        return _delta_t;
    }
    void setDelta_t(double delta_t) {
        BaseDiscretize::_delta_t = delta_t;
    }
    virtual ~BaseDiscretize() = default;

protected:
    double _delta_t;
    normal_distribution<> _normal_gen;
    BaseSDE *_model;
    double random_normal() {
        // this function always sample from a standard normal distribution.
        return _normal_gen(rgen);
    }
};

class EulerDiscretize:BaseDiscretize {
public:
    using BaseDiscretize::BaseDiscretize;

    double next_x(double x) override {
        return x + _model->mu(x) * _delta_t + _model->sigma(x) * sqrt(_delta_t) * random_normal();
    }
};

class MilsteinDiscretize:BaseDiscretize {
public:
    using BaseDiscretize::BaseDiscretize;

    double next_x(double x) override {
        double random_increment = sqrt(_delta_t) * random_normal();

        return x + _model->mu(x) * _delta_t + _model->sigma(x) * random_increment
        + 0.5 * _model->d_sigma_wrt(x) * _model->sigma(x) * (pow(random_increment, 2.0) - _delta_t);
    }
};



#endif //STAT906A3_Q3_H
