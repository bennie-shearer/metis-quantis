// quantlib_credit_risk_model.hpp
#pragma once

#include <ql/quantlib.hpp>
#include <vector>
#include <random>

namespace quant {

class QuantLibCreditRiskModel {
private:
    double defaultCorrelation_;
    size_t numObligors_;
    std::vector<double> defaultProbabilities_;
    std::vector<double> exposures_;
    std::vector<double> recoveryRates_;
    QuantLib::Matrix correlationMatrix_;

    // Private methods
    void initializeCorrelationMatrix();
    std::vector<double> simulateDefaultTimes(
        const QuantLib::Matrix& correlationMatrix,
        const std::vector<double>& hazardRates,
        size_t numSimulations);

public:
    // Constructors
    QuantLibCreditRiskModel();
    QuantLibCreditRiskModel(double defaultCorrelation, size_t numObligors);

    // Configuration methods
    void setDefaultProbabilities(const std::vector<double>& probabilities);
    void setExposures(const std::vector<double>& exposures);
    void setRecoveryRates(const std::vector<double>& recoveryRates);

    // Risk calculation methods
    double calculateExpectedLoss(); // No-parameter version
    double calculateExpectedLoss(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate);

    std::vector<double> simulateLosses(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate,
        size_t numSimulations);

    std::vector<double> simulateDefaultTimes(
        const QuantLib::Date& startDate,
        size_t numSimulations) const;

    double calculateValueAtRisk(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate,
        double confidenceLevel,
        size_t numSimulations);

    double calculateExpectedShortfall(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate,
        double confidenceLevel,
        size_t numSimulations);

    std::vector<double> calculateMarginalVaR(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate,
        double confidenceLevel,
        size_t numSimulations);
};

} // namespace quant