// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_vasicek_model.cpp
#include "quantlib_vasicek_model.hpp"
#include "metis_json.hpp"
#include <ql/quantlib.hpp>
#include <vector>
#include <cmath>
#include <iostream>

using namespace simple_json;

namespace {

// Function to price zero-coupon bond under Vasicek model
double vasicekBondPrice(double r0, double a, double b, double sigma, double T) {
    double B = (1.0 - std::exp(-a * T)) / a;
    double A = std::exp((b - sigma * sigma / (2.0 * a * a)) * (B - T) - sigma * sigma * B * B / (4.0 * a));
    return A * std::exp(-B * r0);
}

// Mean-squared error function for fitting Vasicek model to data
double vasicekMSE(const std::vector<double>& times,
                  const std::vector<double>& rates,
                  double r0, double a, double b, double sigma) {
    double error = 0.0;
    for (size_t i = 0; i < times.size(); ++i) {
        double modelPrice = vasicekBondPrice(r0, a, b, sigma, times[i]);
        double marketRate = rates[i];
        double modelRate = -std::log(modelPrice) / times[i];
        double diff = modelRate - marketRate;
        error += diff * diff;
    }
    return error / times.size();
}

// A simpler cost function for Vasicek model calibration
class VasicekSimpleCostFunction {
private:
    const std::vector<double>& times_;
    const std::vector<double>& rates_;
    double r0_;

public:
    VasicekSimpleCostFunction(const std::vector<double>& times,
                              const std::vector<double>& rates,
                              double r0)
        : times_(times), rates_(rates), r0_(r0) {}

    double operator()(double a, double b, double sigma) const {
        return vasicekMSE(times_, rates_, r0_, a, b, sigma);
    }
};

// Numerical optimization to find Vasicek parameters
std::tuple<double, double, double> optimizeVasicekParameters(
    const std::vector<double>& times,
    const std::vector<double>& rates,
    double r0) {

    // Initial guess for parameters
    double a_best = 0.1;
    double b_best = 0.05;
    double sigma_best = 0.01;
    double best_error = vasicekMSE(times, rates, r0, a_best, b_best, sigma_best);

    // Simple grid search for optimization
    for (double a = 0.01; a <= 0.5; a += 0.01) {
        for (double b = 0.01; b <= 0.1; b += 0.005) {
            for (double sigma = 0.001; sigma <= 0.05; sigma += 0.001) {
                double error = vasicekMSE(times, rates, r0, a, b, sigma);
                if (error < best_error) {
                    a_best = a;
                    b_best = b;
                    sigma_best = sigma;
                    best_error = error;
                }
            }
        }
    }

    return {a_best, b_best, sigma_best};
}

} // anonymous namespace

namespace quant {

SimpleJSON calibrateVasicekModel(const SimpleJSON& params) {
    SimpleJSON result = SimpleJSON::object();

    try {
        // Extract parameters
        const auto& timeArray = params["times"].asArray();
        const auto& rateArray = params["rates"].asArray();
        double r0 = params["r0"].asDouble();

        // Convert JSON arrays to vectors
        std::vector<double> times;
        std::vector<double> rates;

        for (const auto& time : timeArray) {
            times.push_back(time.asDouble());
        }

        for (const auto& rate : rateArray) {
            rates.push_back(rate.asDouble());
        }

        // Calibrate model parameters
        auto [a, b, sigma] = optimizeVasicekParameters(times, rates, r0);

        // Calculate error
        double error = vasicekMSE(times, rates, r0, a, b, sigma);

        // Populate result
        result["a"] = a;
        result["b"] = b;
        result["sigma"] = sigma;
        result["r0"] = r0;
        result["error"] = error;

        // Return model-implied rates
        SimpleJSON modelRates = SimpleJSON::array();
        for (size_t i = 0; i < times.size(); ++i) {
            double modelPrice = vasicekBondPrice(r0, a, b, sigma, times[i]);
            double modelRate = -std::log(modelPrice) / times[i];
            modelRates.pushBack(modelRate);
        }
        result["modelRates"] = modelRates;

        // Create a new array for times in the response
        SimpleJSON timeArrayOut = SimpleJSON::array();
        for (size_t i = 0; i < times.size(); ++i) {
            timeArrayOut.pushBack(times[i]);
        }

        // Use the array for output
        result["times"] = timeArrayOut;

    } catch (const std::exception& e) {
        result["error"] = e.what();
    }

    return result;
}

SimpleJSON simulateVasicekRates(const SimpleJSON& params) {
    SimpleJSON result = SimpleJSON::object();

    try {
        // Extract parameters
        double a = params["a"].asDouble();
        double b = params["b"].asDouble();
        double sigma = params["sigma"].asDouble();
        double r0 = params["r0"].asDouble();
        double T = params["horizon"].asDouble();
        int steps = params["steps"].asInteger();
        int paths = params["paths"].asInteger();

        // Time step
        double dt = T / steps;

        // Random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> normal(0.0, 1.0);

        // Simulate paths
        std::vector<std::vector<double>> simulations(paths, std::vector<double>(steps + 1));

        for (int i = 0; i < paths; ++i) {
            simulations[i][0] = r0;

            for (int j = 0; j < steps; ++j) {
                double r = simulations[i][j];
                double drift = a * (b - r) * dt;
                double diffusion = sigma * std::sqrt(dt) * normal(gen);
                simulations[i][j + 1] = r + drift + diffusion;
            }
        }

        // Convert to JSON format
        SimpleJSON pathsArray = SimpleJSON::array();

        for (int i = 0; i < paths; ++i) {
            SimpleJSON path = SimpleJSON::array();
            for (int j = 0; j <= steps; ++j) {
                path.pushBack(simulations[i][j]);
            }
            pathsArray.pushBack(path);
        }

        // Create time points
        SimpleJSON timePoints = SimpleJSON::array();
        for (int j = 0; j <= steps; ++j) {
            timePoints.pushBack(j * dt);
        }

        // Populate result
        result["paths"] = pathsArray;
        result["times"] = timePoints;
        result["parameters"] = SimpleJSON::object();
        result["parameters"]["a"] = a;
        result["parameters"]["b"] = b;
        result["parameters"]["sigma"] = sigma;
        result["parameters"]["r0"] = r0;

    } catch (const std::exception& e) {
        result["error"] = e.what();
    }

    return result;
}

// Calculate ZCB price using the Vasicek model
double calculateVasicekBondPrice(double r0, double a, double b, double sigma, double T) {
    return vasicekBondPrice(r0, a, b, sigma, T);
}

} // namespace quant