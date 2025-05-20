// Path: include/quantlib_heston_model.hpp

#ifndef QUANTLIB_HESTON_MODEL_H
#define QUANTLIB_HESTON_MODEL_H

#include <ql/quantlib.hpp>
#include <vector>

namespace quant {

class QuantLibHestonModel {
public:
    struct HestonParameters {
        double v0;       // Initial variance
        double kappa;    // Mean reversion speed of variance
        double theta;    // Long-term variance
        double sigma;    // Volatility of variance
        double rho;      // Correlation between asset and variance
    };

    struct CalibrationResults {
        HestonParameters params;
        double objectiveValue;
        bool success;
        std::vector<double> modelPrices;
        std::vector<double> marketPrices;
        std::vector<double> relativeErrors;
    };

    // Constructor with default parameters
    QuantLibHestonModel(double v0 = 0.04,
                       double kappa = 1.0,
                       double theta = 0.04,
                       double sigma = 0.2,
                       double rho = -0.5)
        : m_params({v0, kappa, theta, sigma, rho}) {}

    // Create Heston process
    QuantLib::ext::shared_ptr<QuantLib::HestonProcess> createHestonProcess(
        double spotPrice,
        double riskFreeRate,
        double dividendYield,
        double timeToMaturity);

    // Create Heston model
    QuantLib::ext::shared_ptr<QuantLib::HestonModel> createHestonModel(
        const QuantLib::ext::shared_ptr<QuantLib::HestonProcess>& process);

    // Price European option with Heston model
    double priceEuropeanOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call);

    // Price European options for multiple strikes and maturities
    std::vector<double> priceEuropeanOptionGrid(
        double spotPrice,
        const std::vector<double>& strikes,
        const std::vector<double>& maturities,
        double riskFreeRate,
        double dividendYield,
        QuantLib::Option::Type optionType = QuantLib::Option::Call);

    // Calculate implied volatility surface using Heston model
    std::vector<std::vector<double>> calculateImpliedVolatilitySurface(
        double spotPrice,
        const std::vector<double>& strikes,
        const std::vector<double>& maturities,
        double riskFreeRate,
        double dividendYield);

    // Calibrate Heston model to market data
    CalibrationResults calibrateHestonModel(
        double spotPrice,
        const std::vector<double>& strikes,
        const std::vector<double>& maturities,
        const std::vector<double>& marketVolatilities,
        double riskFreeRate,
        double dividendYield);

    // Get/set Heston parameters
    HestonParameters getParameters() const { return m_params; }
    void setParameters(const HestonParameters& params) { m_params = params; }

private:
    HestonParameters m_params;
};

} // namespace quant

#endif // QUANTLIB_HESTON_MODEL_H