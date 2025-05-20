// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_risk_management.hpp
#ifndef QUANTLIB_RISK_MANAGEMENT_H
#define QUANTLIB_RISK_MANAGEMENT_H

#include <ql/quantlib.hpp>
#include <vector>
#include <string>

namespace quant {

class QuantLibRiskManagement {
public:
    // Results structure for VaR calculations
    struct VaRResult {
        double varPercent = 0.0;      // VaR as percentage of portfolio
        double varAmount = 0.0;       // VaR in monetary terms
        double esPercent = 0.0;       // Expected Shortfall (CVaR) as percentage
        double esAmount = 0.0;        // Expected Shortfall in monetary terms
    };

    // Stress scenario structure
    struct StressScenario {
        std::string name;
        std::vector<double> assetImpacts;  // Percentage impacts on each asset
    };

    // Stress test result structure
    struct StressTestResult {
        std::string scenarioName;
        double portfolioPercentChange = 0.0;
        double portfolioValueChange = 0.0;
        double newPortfolioValue = 0.0;
    };

    // Historical simulation VaR
    VaRResult calculateVaR(
        const std::vector<double>& historicalReturns,
        double confidenceLevel = 0.95,
        double portfolioValue = 1000000.0);

    // Parametric VaR
    VaRResult calculateParametricVaR(
        double portfolioReturn,
        double portfolioVolatility,
        double confidenceLevel = 0.95,
        double portfolioValue = 1000000.0,
        int timeHorizon = 1);  // Days

    // Monte Carlo VaR
    VaRResult calculateMonteCarloVaR(
        const std::vector<double>& expectedReturns,
        const QuantLib::Matrix& covarianceMatrix,
        const std::vector<double>& weights,
        double confidenceLevel = 0.95,
        double portfolioValue = 1000000.0,
        int timeHorizon = 1,  // Days
        size_t numSimulations = 10000);

    // Stress testing
    std::vector<StressTestResult> performStressTest(
        const std::vector<StressScenario>& scenarios,
        const std::vector<double>& weights,
        double portfolioValue = 1000000.0);

    // Portfolio volatility calculation
    double calculatePortfolioVolatility(
        const QuantLib::Matrix& covarianceMatrix,
        const std::vector<double>& weights);
};

} // namespace quant

#endif // QUANTLIB_RISK_MANAGEMENT_H