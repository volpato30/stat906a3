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

class ShortInterestModel:public BaseSDE {
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
        return _std * _gamma * pow(x, _gamma - 1.);
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

class EulerDiscretize: public BaseDiscretize {
public:
    using BaseDiscretize::BaseDiscretize;

    double next_x(double x) override {
        return x + _model->mu(x) * _delta_t + _model->sigma(x) * sqrt(_delta_t) * random_normal();
    }
};

class MilsteinDiscretize: public BaseDiscretize {
public:
    using BaseDiscretize::BaseDiscretize;

    double next_x(double x) override {
        double random_increment = sqrt(_delta_t) * random_normal();

        return x + _model->mu(x) * _delta_t + _model->sigma(x) * random_increment
        + 0.5 * _model->d_sigma_wrt(x) * _model->sigma(x) * (pow(random_increment, 2.0) - _delta_t);
    }
};


template <class T>
T mean(const std::vector<T> x) {
    T result=0;
    for (const T &c : x) {
        result += c;
    }
    result = result / x.size();
    return result;
}

template <class T>
T standard_deviation(const std::vector<T> x) {
    T result=0;
    T m = mean(x);
    for (const T &c : x) {
        result += (c - m) * (c - m);
    }
    return sqrt(result / (x.size() - 1));
}


template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
    return (1 - t)*v0 + t*v1;
}

double empirical_lower_quantile(std::vector<double> &x, double lower_quantile) {
    // inspired by https://stackoverflow.com/questions/11964552/finding-quartiles
    double position = Lerp<double>(-0.5, x.size() - 0.5, lower_quantile);
    auto left = static_cast<unsigned long>(std::max(floor(position), 0.));
    auto right = std::min(static_cast<unsigned long>(ceil(position)), x.size());

    std::nth_element (x.begin(), x.begin()+left, x.end());
    double low_ = x.at(left);

    std::nth_element (x.begin(), x.begin()+right, x.end());
    double high_ = x.at(right);

    return Lerp<double>(low_, high_, position - left);

}

double exp(BaseDiscretize *disc, double r0, int n_path=100) {
    std::vector<double> x_paths;
    for (int i=0; i < n_path; i++) {
        x_paths.push_back(disc->get_x_t(6.0, r0));
    }
    return empirical_lower_quantile(x_paths, 0.95);
}

void q3() {
    double r0 = 0.06, b = 0.0654, kappa = 0.1863, std = 0.08, gamma = 0.5;
    int n_exp = 100;
    ShortInterestModel cirModel(kappa, b, std, gamma);
    EulerDiscretize eular(1/250., &cirModel);
    MilsteinDiscretize milstein(1/250., &cirModel);

    std::cout << "delta t: " << 1 / 250. << std::endl;

    std::vector<double> result;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < n_exp; i++) {
        result.push_back(exp(&eular, r0));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Eular:\t" << "mean: " << mean(result) << "\tstd: " << standard_deviation(result) <<
    "\t takes " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "miliseconds" << std::endl;

    result.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < n_exp; i++) {
        result.push_back(exp(&milstein, r0));
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "milstein:\t" << "mean: " << mean(result) << "\tstd: " << standard_deviation(result) <<
    "\t takes " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "miliseconds" << std::endl;

    std::cout << "delta t: " << 1 / 4. << std::endl;
    eular.setDelta_t(1/4.);
    milstein.setDelta_t(1/4.);

    result.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < n_exp; i++) {
        result.push_back(exp(&eular, r0));
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Eular:\t" << "mean: " << mean(result) << "\tstd: " << standard_deviation(result) <<
              "\t takes " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " miliseconds" << std::endl;

    result.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < n_exp; i++) {
        result.push_back(exp(&milstein, r0));
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "milstein:\t" << "mean: " << mean(result) << "\tstd: " << standard_deviation(result) <<
              "\t takes " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " miliseconds" << std::endl;



}

#endif //STAT906A3_Q3_H
