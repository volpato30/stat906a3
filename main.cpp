//
// Created by rui on 11/12/18.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::random_device rd;
std::mt19937 rgen(rd());
using std::normal_distribution;
using std::uniform_real_distribution;


struct MCMCResult {
    std::vector<double> x1Chain;
    std::vector<double> x2Chain;
    double rejection_rate=0;
};

double gaussian_pdf(double &x, double mu, double sigma) {
    return ( 1 / ( sigma * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (x-mu)/sigma, 2.0 ) );
}

double q_pdf(double &x1, double &x2, double &sigma) {
    /// compute the proposal densitiy
    return gaussian_pdf(x1, 0, sigma) * gaussian_pdf(x1, 0, sigma);
}

double target_pdf(double &x1, double &x2, double &theta) {
    /// Use clayton copula as target distribution
    if (x1 > 1 || x1 < 0) {
        return 0;
    }
    if (x2 > 1 || x2 < 0) {
        return 0;
    }
    return (theta + 1) * pow(x1 * x2, -theta-1) * pow(pow(x1, -theta) + pow(x2, -theta) - 1, -(2*theta+1)/theta);
}


MCMCResult MCMC(unsigned long num_samples, double sigma) {
    double theta = 8;
    normal_distribution<> normal_gen(0, sigma);
    uniform_real_distribution<double> unif_gen(0, 1);
    double x1_proposal, x2_proposal, acceptance_rate, trial, total_rejection_rate=0;

    // given the property of target distribution, suppose we start from point (0.5, 0.5)
    // generate proposal from independent gaussian distribution.
    MCMCResult result={};
    result.x1Chain.clear();
    result.x2Chain.clear();
    result.x1Chain.push_back(0.5);
    result.x2Chain.push_back(0.5);

    for (unsigned long i=0; i<num_samples; i++) {
        x1_proposal = result.x1Chain.at(i) + normal_gen(rgen);
        x2_proposal = result.x2Chain.at(i) + normal_gen(rgen);
        acceptance_rate = std::min(1.0, target_pdf(x1_proposal, x2_proposal, theta) /
                                        target_pdf(result.x1Chain.at(i), result.x2Chain.at(i), theta));

        trial = unif_gen(rgen);
        if (trial < acceptance_rate) {
            result.x1Chain.push_back(x1_proposal);
            result.x2Chain.push_back(x2_proposal);
        } else {
            result.x1Chain.push_back(result.x1Chain.at(i));
            result.x2Chain.push_back(result.x2Chain.at(i));
            total_rejection_rate += 1;
        }
    }
    result.rejection_rate = total_rejection_rate / double(num_samples);
    return result;
}

void q2(){
//    std::vector<double> sigmaList = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
//    MCMCResult mc_result;
//    for (const double &sigma : sigmaList) {
//        mc_result = MCMC(100000, sigma);
//        std::cout << "sigma: " << sigma << "\trejection rate: " << mc_result.rejection_rate << std::endl;
//    }
    // thus we choose sigma=0.2 because it has an acceptance rate of around 20%
    unsigned long num_samples = 30000, i=0;
    std::vector<double> index(num_samples + 1);
    for (i=0; i < num_samples+1; i++) {
        index.at(i) = double(i);
    }
    std::vector<std::string> labels = {"path 1", "path 2", "path 3", "path 4", "path 5"};
    std::map<std::string, std::string> keyword;


    std::vector<MCMCResult> mc_results(5);
    for (i = 0; i<5; i++) {
        mc_results.at(i) = MCMC(num_samples, 0.2);
    }
    // draw path
//    plt::figure_size(200, 500);
    for (i=0; i<5; i++) {
        keyword["label"] = labels.at(i);
        plt::subplot(5, 1, i+1);
        plt::plot(index, mc_results.at(i).x1Chain);
        plt::title(labels.at(i));
    }
    plt::save("x1_5_path.png");

    plt::figure();
    for (i=0; i<5; i++) {
        keyword["label"] = labels.at(i);
        plt::subplot(5, 1, i+1);
        plt::plot(index, mc_results.at(i).x2Chain);
        plt::title(labels.at(i));
    }
    plt::save("x2_5_path.png");

    unsigned long warm_up = 10000;
    std::vector<double> x1(&mc_results.at(0).x1Chain[warm_up], &mc_results.at(0).x1Chain[num_samples]);
    std::vector<double> x2(&mc_results.at(0).x2Chain[warm_up], &mc_results.at(0).x2Chain[num_samples]);

    plt::figure();
    plt::plot(x1, x2, "r.");
    plt::title("MCMC generated sample points for clayton copula where theta=8");
    plt::save("scatter.png");

}


int main() {
    q2();
}