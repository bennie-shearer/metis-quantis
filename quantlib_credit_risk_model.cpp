#include <random>
#include "quantlib_credit_risk_model.hpp"

using namespace QuantLib;

namespace quant {

QuantLibCreditRiskModel::QuantLibCreditRiskModel()
    : defaultCorrelation_(0.0), numObligors_(0) {
}

QuantLibCreditRiskModel::QuantLibCreditRiskModel(
    double defaultCorrelation,
    size_t numObligors)
    : defaultCorrelation_(defaultCorrelation),
      numObligors_(numObligors) {

    // Initialize vectors with the right size
    defaultProbabilities_.resize(numObligors_, 0.0);
    exposures_.resize(numObligors_, 0.0);
    recoveryRates_.resize(numObligors_, 0.4); // Default recovery rate

    // Initialize correlation matrix
    initializeCorrelationMatrix();
}

void QuantLibCreditRiskModel::initializeCorrelationMatrix() {
    // Create the correlation matrix with the specified default correlation
    correlationMatrix_ = Matrix(numObligors_, numObligors_, 0.0);

    // Set diagonal to 1.0 and off-diagonal to defaultCorrelation_
    for (Size i = 0; i < numObligors_; ++i) {
        correlationMatrix_[i][i] = 1.0;
        for (Size j = 0; j < i; ++j) {
            correlationMatrix_[i][j] = correlationMatrix_[j][i] = defaultCorrelation_;
        }
    }
}

void QuantLibCreditRiskModel::setDefaultProbabilities(const std::vector<double>& probabilities) {
    if (probabilities.size() != numObligors_) {
        // Handle error: wrong size
        throw std::invalid_argument("Size of default probabilities must match number of obligors");
    }
    defaultProbabilities_ = probabilities;
}

void QuantLibCreditRiskModel::setExposures(const std::vector<double>& exposures) {
    if (exposures.size() != numObligors_) {
        // Handle error: wrong size
        throw std::invalid_argument("Size of exposures must match number of obligors");
    }
    exposures_ = exposures;
}

void QuantLibCreditRiskModel::setRecoveryRates(const std::vector<double>& recoveryRates) {
    if (recoveryRates.size() != numObligors_) {
        // Handle error: wrong size
        throw std::invalid_argument("Size of recovery rates must match number of obligors");
    }
    recoveryRates_ = recoveryRates;
}

double QuantLibCreditRiskModel::calculateExpectedLoss() {
    double expectedLoss = 0.0;

    // Calculate the expected loss for each obligor
    for (Size i = 0; i < numObligors_; ++i) {
        double pd = defaultProbabilities_[i];                  // Probability of default
        double exposure = exposures_[i];                       // Exposure at default
        double lgd = 1.0 - recoveryRates_[i];                  // Loss given default

        expectedLoss += pd * exposure * lgd;
    }

    return expectedLoss;
}

double QuantLibCreditRiskModel::calculateExpectedLoss(
    const QuantLib::Date& startDate,
    const QuantLib::Date& endDate) {

    // Calculate expected loss between dates
    // This method could use the term structure to adjust probabilities
    // For now, just delegate to the simpler method
    return calculateExpectedLoss();
}

std::vector<double> QuantLibCreditRiskModel::simulateDefaultTimes(
    const QuantLib::Matrix& correlationMatrix,
    const std::vector<double>& hazardRates,
    size_t numSimulations) {

    // Create a random number generator
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::normal_distribution<> stdNormal(0.0, 1.0);

    // Simulate default times
    std::vector<double> defaultTimes(numObligors_ * numSimulations);

    // Perform Cholesky decomposition of the correlation matrix
    Matrix cholMatrix = CholeskyDecomposition(correlationMatrix);

    // Simulate default times
    for (Size sim = 0; sim < numSimulations; ++sim) {
        // Generate standard normal variables
        std::vector<double> z(numObligors_);
        for (Size i = 0; i < numObligors_; ++i) {
            z[i] = stdNormal(gen);
        }

        // Apply Cholesky matrix to get correlated normals
        std::vector<double> correlated(numObligors_, 0.0);
        for (Size i = 0; i < numObligors_; ++i) {
            for (Size j = 0; j <= i; ++j) {
                correlated[i] += cholMatrix[i][j] * z[j];
            }
        }

        // Convert to default times
        for (Size i = 0; i < numObligors_; ++i) {
            // Calculate default time using inverse CDF of exponential distribution
            double u = CumulativeNormalDistribution()(correlated[i]);
            defaultTimes[sim * numObligors_ + i] = -std::log(1.0 - u) / hazardRates[i];
        }
    }

    return defaultTimes;
}

std::vector<double> QuantLibCreditRiskModel::simulateDefaultTimes(
    const QuantLib::Date& startDate,
    size_t numSimulations) const {

    // Implement the const version required by the header
    // You might want to adjust this implementation as needed

    // Convert default probabilities to hazard rates
    std::vector<double> hazardRates(numObligors_);
    for (Size i = 0; i < numObligors_; ++i) {
        double pd = defaultProbabilities_[i];
        // Convert PD to hazard rate (simple approximation)
        hazardRates[i] = -std::log(1.0 - pd);
    }

    // Call the non-const version with correlationMatrix_
    // Note: This is a workaround for demonstration purposes
    // In real code, you'd implement this independently
    return const_cast<QuantLibCreditRiskModel*>(this)->simulateDefaultTimes(
        correlationMatrix_, hazardRates, numSimulations);
}

std::vector<double> QuantLibCreditRiskModel::simulateLosses(
    const QuantLib::Date& startDate,
    const QuantLib::Date& endDate,
    size_t numSimulations) {

    // Convert default probabilities to hazard rates
    std::vector<double> hazardRates(numObligors_);
    for (Size i = 0; i < numObligors_; ++i) {
        double pd = defaultProbabilities_[i];
        // Convert PD to hazard rate (simple approximation)
        hazardRates[i] = -std::log(1.0 - pd);
    }

    // Simulate default times
    std::vector<double> defaultTimes = simulateDefaultTimes(correlationMatrix_, hazardRates, numSimulations);

    // Calculate horizon in years (approximate)
    double horizonYears = Real(endDate - startDate) / 365.0;

    // Calculate losses for each simulation
    std::vector<double> losses(numSimulations, 0.0);

    for (Size sim = 0; sim < numSimulations; ++sim) {
        double portfolioLoss = 0.0;

        for (Size i = 0; i < numObligors_; ++i) {
            double defaultTime = defaultTimes[sim * numObligors_ + i];

            // Check if default occurs within the horizon
            if (defaultTime <= horizonYears) {
                double exposureAtDefault = exposures_[i];
                double lossGivenDefault = exposureAtDefault * (1.0 - recoveryRates_[i]);
                portfolioLoss += lossGivenDefault;
            }
        }

        losses[sim] = portfolioLoss;
    }

    return losses;
}

double QuantLibCreditRiskModel::calculateValueAtRisk(
    const QuantLib::Date& startDate,
    const QuantLib::Date& endDate,
    double confidenceLevel,
    size_t numSimulations) {

    // Simulate losses
    std::vector<double> losses = simulateLosses(startDate, endDate, numSimulations);

    // Sort losses
    std::sort(losses.begin(), losses.end());

    // Calculate VaR
    size_t index = static_cast<size_t>(std::ceil(numSimulations * (1.0 - confidenceLevel)) - 1);
    if (index >= losses.size()) {
        index = losses.size() - 1;
    }

    return losses[index];
}

double QuantLibCreditRiskModel::calculateExpectedShortfall(
    const QuantLib::Date& startDate,
    const QuantLib::Date& endDate,
    double confidenceLevel,
    size_t numSimulations) {

    // Simulate losses
    std::vector<double> losses = simulateLosses(startDate, endDate, numSimulations);

    // Sort losses
    std::sort(losses.begin(), losses.end());

    // Calculate VaR index
    size_t VaRIndex = static_cast<size_t>(std::ceil(numSimulations * (1.0 - confidenceLevel)) - 1);
    if (VaRIndex >= losses.size()) {
        VaRIndex = losses.size() - 1;
    }

    // Calculate Expected Shortfall (average of losses beyond VaR)
    double ES = 0.0;
    for (size_t i = VaRIndex; i < losses.size(); ++i) {
        ES += losses[i];
    }

    ES /= (losses.size() - VaRIndex);

    return ES;
}

std::vector<double> QuantLibCreditRiskModel::calculateMarginalVaR(
    const QuantLib::Date& startDate,
    const QuantLib::Date& endDate,
    double confidenceLevel,
    size_t numSimulations) {

    // Calculate portfolio VaR
    double portfolioVaR = calculateValueAtRisk(startDate, endDate, confidenceLevel, numSimulations);

    // Calculate marginal VaR for each obligor
    std::vector<double> marginalVaRs(numObligors_, 0.0);

    // Use a small increment for numerical derivative
    double epsilon = 0.001;

    for (Size i = 0; i < numObligors_; ++i) {
        // Save original exposure
        double originalExposure = exposures_[i];

        // Increase exposure by a small amount
        exposures_[i] *= (1.0 + epsilon);

        // Calculate new VaR
        double newVaR = calculateValueAtRisk(startDate, endDate, confidenceLevel, numSimulations);

        // Calculate marginal VaR (numerical derivative)
        marginalVaRs[i] = (newVaR - portfolioVaR) / (originalExposure * epsilon);

        // Restore original exposure
        exposures_[i] = originalExposure;
    }

    return marginalVaRs;
}

} // namespace quant