// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_var_calculator.cpp
#include <random>
#include "quantlib_var_calculator.hpp"
#include "metis_json.hpp"
#include <ql/math/statistics/sequencestatistics.hpp>
#include <ql/math/statistics/generalstatistics.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/matrix.hpp>

using namespace simple_json;

namespace quant {

double calculateHistoricalVaR(const std::vector<double>& returns, double confidence) {
    if (returns.empty()) {
        return 0.0;
    }

    // Copy returns for sorting
    std::vector<double> sortedReturns = returns;
    std::sort(sortedReturns.begin(), sortedReturns.end());

    // Calculate VaR index
    int varIndex = static_cast<int>(sortedReturns.size() * (1.0 - confidence));

    // Ensure index is within bounds
    varIndex = std::max(0, std::min(varIndex, static_cast<int>(sortedReturns.size() - 1)));

    // Historical VaR is the negative of the return at the VaR index
    return -sortedReturns[varIndex];
}

double calculateParametricVaR(double mean, double stddev, double confidence) {
    // Create normal distribution
    QuantLib::InverseCumulativeNormal invCumNormal;

    // Calculate VaR
    double z = invCumNormal(1.0 - confidence);
    return -(mean + z * stddev);
}

double calculateMonteCarloVaR(
    const std::vector<double>& returns,
    double confidence,
    int numSimulations,
    int horizon) {

    if (returns.empty() || numSimulations <= 0 || horizon <= 0) {
        return 0.0;
    }

    // Calculate mean and standard deviation of returns
    double sum = 0.0;
    for (double r : returns) {
        sum += r;
    }
    double mean = sum / returns.size();

    double sumSq = 0.0;
    for (double r : returns) {
        sumSq += (r - mean) * (r - mean);
    }
    double variance = sumSq / (returns.size() - 1);
    double stddev = std::sqrt(variance);

    // Set up random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal(mean, stddev);

    // Simulate returns
    std::vector<double> simulatedReturns;
    simulatedReturns.reserve(numSimulations);

    for (int i = 0; i < numSimulations; ++i) {
        double simulatedReturn = 0.0;

        // Simulate return over the horizon
        for (int j = 0; j < horizon; ++j) {
            simulatedReturn += normal(gen);
        }

        simulatedReturns.push_back(simulatedReturn);
    }

    // Calculate VaR from simulated returns
    return calculateHistoricalVaR(simulatedReturns, confidence);
}

SimpleJSON calculateVaR(const SimpleJSON& params) {
    SimpleJSON result = SimpleJSON::object();

    try {
        // Extract parameters
        double confidence = params["confidence"].asDouble();
        std::string method = params["method"].asString();

        // Extract returns
        std::vector<double> returns;
        for (const auto& r : params["returns"].asArray()) {
            returns.push_back(r.asDouble());
        }

        // Calculate VaR using the specified method
        double var = 0.0;

        if (method == "historical") {
            var = calculateHistoricalVaR(returns, confidence);
        } else if (method == "parametric") {
            // Calculate mean and standard deviation
            double sum = 0.0;
            for (double r : returns) {
                sum += r;
            }
            double mean = sum / returns.size();

            double sumSq = 0.0;
            for (double r : returns) {
                sumSq += (r - mean) * (r - mean);
            }
            double variance = sumSq / (returns.size() - 1);
            double stddev = std::sqrt(variance);

            var = calculateParametricVaR(mean, stddev, confidence);
        } else if (method == "montecarlo") {
            int numSimulations = params["simulations"].asInteger();
            int timeHorizon = params["timeHorizon"].asInteger();
            var = calculateMonteCarloVaR(returns, confidence, numSimulations, timeHorizon);
        } else {
            throw std::runtime_error("Unknown VaR calculation method");
        }

        // Populate result
        result["var"] = var;
        result["confidence"] = confidence;
        result["method"] = method;

        if (method == "montecarlo") {
            result["simulations"] = params["simulations"].asInteger();
            // Fix for error: Create an integer assignment operator for SimpleJSON
            double timeHorizon = params["timeHorizon"].asDouble();
            result["timeHorizon"] = timeHorizon;
        }
    } catch (const std::exception& e) {
        result["error"] = e.what();
    }

    return result;
}

} // namespace quant