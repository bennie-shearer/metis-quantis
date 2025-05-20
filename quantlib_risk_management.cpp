// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_risk_management.cpp
#include "quantlib_risk_management.hpp"
#include <ql/math/statistics/generalstatistics.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/math/randomnumbers/sobolrsg.hpp>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace quant {

QuantLibRiskManagement::VaRResult QuantLibRiskManagement::calculateVaR(
    const std::vector<double>& historicalReturns,
    double confidenceLevel,
    double portfolioValue) {

    try {
        VaRResult result;

        // Sort returns from worst to best (lowest to highest)
        std::vector<double> sortedReturns = historicalReturns;
        std::sort(sortedReturns.begin(), sortedReturns.end());

        // Calculate VaR index based on confidence level
        size_t varIndex = static_cast<size_t>(std::floor(sortedReturns.size() * (1.0 - confidenceLevel)));
        if (varIndex >= sortedReturns.size()) {
            varIndex = sortedReturns.size() - 1;
        }

        // VaR is the loss at the specified confidence level
        result.varPercent = -sortedReturns[varIndex];
        result.varAmount = portfolioValue * result.varPercent;

        // Calculate Expected Shortfall (Conditional VaR)
        double esSum = 0.0;
        for (size_t i = 0; i < varIndex; ++i) {
            esSum += -sortedReturns[i];
        }

        result.esPercent = (varIndex > 0) ? esSum / varIndex : result.varPercent;
        result.esAmount = portfolioValue * result.esPercent;

        return result;
    } catch (std::exception& e) {
        std::cerr << "Error calculating VaR: " << e.what() << std::endl;
        return VaRResult();
    }
}

QuantLibRiskManagement::VaRResult QuantLibRiskManagement::calculateParametricVaR(
    double portfolioReturn,
    double portfolioVolatility,
    double confidenceLevel,
    double portfolioValue,
    int timeHorizon) {

    try {
        VaRResult result;

        // Calculate confidence factor using inverse normal distribution
        QuantLib::InverseCumulativeNormal inverseNormal;
        double alpha = inverseNormal(confidenceLevel);

        // Calculate VaR for the specified time horizon
        double varPercent = -portfolioReturn + alpha * portfolioVolatility * std::sqrt(static_cast<double>(timeHorizon) / 252.0);
        result.varPercent = varPercent;
        result.varAmount = portfolioValue * varPercent;

        // For normal distribution, ES = VaR + (pdf(alpha) / (1-CL)) * sigma
        QuantLib::CumulativeNormalDistribution normal;
        double pdf = normal.derivative(alpha);
        result.esPercent = varPercent + (pdf / (1.0 - confidenceLevel)) * portfolioVolatility * std::sqrt(static_cast<double>(timeHorizon) / 252.0);
        result.esAmount = portfolioValue * result.esPercent;

        return result;
    } catch (std::exception& e) {
        std::cerr << "Error calculating parametric VaR: " << e.what() << std::endl;
        return VaRResult();
    }
}

QuantLibRiskManagement::VaRResult QuantLibRiskManagement::calculateMonteCarloVaR(
    const std::vector<double>& expectedReturns,
    const QuantLib::Matrix& covarianceMatrix,
    const std::vector<double>& weights,
    double confidenceLevel,
    double portfolioValue,
    int timeHorizon,
    size_t numSimulations) {

    try {
        VaRResult result;

        // Check for valid input dimensions
        if (expectedReturns.size() != weights.size() ||
            covarianceMatrix.rows() != expectedReturns.size() ||
            covarianceMatrix.columns() != expectedReturns.size()) {
            throw std::invalid_argument("Inconsistent input dimensions");
        }

        // Calculate Cholesky decomposition of covariance matrix
        QuantLib::Matrix cholMatrix = QuantLib::CholeskyDecomposition(covarianceMatrix);

        // Set up random number generator
        QuantLib::MersenneTwisterUniformRng uniformGenerator(42);
        QuantLib::InverseCumulativeNormal icnd;

        // Generate simulated returns
        std::vector<double> simulatedReturns;
        simulatedReturns.reserve(numSimulations);

        for (size_t sim = 0; sim < numSimulations; ++sim) {
            // Generate correlated normal random variables
            std::vector<double> randomFactors(weights.size());
            for (size_t i = 0; i < weights.size(); ++i) {
                randomFactors[i] = icnd(uniformGenerator.next().value);
            }

            // Apply Cholesky to get correlated random samples
            std::vector<double> correlatedReturns(weights.size(), 0.0);
            for (size_t i = 0; i < weights.size(); ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    correlatedReturns[i] += cholMatrix[i][j] * randomFactors[j];
                }
            }

            // Calculate portfolio return for this simulation
            double portfolioReturn = 0.0;
            for (size_t i = 0; i < weights.size(); ++i) {
                double assetReturn = expectedReturns[i] + correlatedReturns[i] *
                                    std::sqrt(static_cast<double>(timeHorizon) / 252.0);
                portfolioReturn += weights[i] * assetReturn;
            }

            simulatedReturns.push_back(portfolioReturn);
        }

        // Sort returns to calculate VaR
        std::sort(simulatedReturns.begin(), simulatedReturns.end());

        // Calculate VaR at specified confidence level
        size_t varIndex = static_cast<size_t>(std::floor(numSimulations * (1.0 - confidenceLevel)));
        if (varIndex >= numSimulations) {
            varIndex = numSimulations - 1;
        }

        result.varPercent = -simulatedReturns[varIndex];
        result.varAmount = portfolioValue * result.varPercent;

        // Calculate Expected Shortfall
        double esSum = 0.0;
        for (size_t i = 0; i < varIndex; ++i) {
            esSum += -simulatedReturns[i];
        }

        result.esPercent = (varIndex > 0) ? esSum / varIndex : result.varPercent;
        result.esAmount = portfolioValue * result.esPercent;

        return result;
    } catch (std::exception& e) {
        std::cerr << "Error calculating Monte Carlo VaR: " << e.what() << std::endl;
        return VaRResult();
    }
}

std::vector<QuantLibRiskManagement::StressTestResult> QuantLibRiskManagement::performStressTest(
    const std::vector<QuantLibRiskManagement::StressScenario>& scenarios,
    const std::vector<double>& weights,
    double portfolioValue) {

    try {
        std::vector<StressTestResult> results;

        for (const auto& scenario : scenarios) {
            StressTestResult result;
            result.scenarioName = scenario.name;

            // Calculate the portfolio impact
            double portfolioImpact = 0.0;
            for (size_t i = 0; i < weights.size() && i < scenario.assetImpacts.size(); ++i) {
                portfolioImpact += weights[i] * scenario.assetImpacts[i];
            }

            result.portfolioPercentChange = portfolioImpact;
            result.portfolioValueChange = portfolioValue * portfolioImpact;
            result.newPortfolioValue = portfolioValue * (1.0 + portfolioImpact);

            results.push_back(result);
        }

        return results;
    } catch (std::exception& e) {
        std::cerr << "Error performing stress test: " << e.what() << std::endl;
        return std::vector<StressTestResult>();
    }
}

double QuantLibRiskManagement::calculatePortfolioVolatility(
    const QuantLib::Matrix& covarianceMatrix,
    const std::vector<double>& weights) {

    try {
        // Calculate portfolio variance
        double portfolioVariance = 0.0;

        for (size_t i = 0; i < weights.size(); ++i) {
            for (size_t j = 0; j < weights.size(); ++j) {
                portfolioVariance += weights[i] * weights[j] * covarianceMatrix[i][j];
            }
        }

        // Return portfolio volatility (standard deviation)
        return std::sqrt(portfolioVariance);
    } catch (std::exception& e) {
        std::cerr << "Error calculating portfolio volatility: " << e.what() << std::endl;
        return -1.0;
    }
}

} // namespace quant