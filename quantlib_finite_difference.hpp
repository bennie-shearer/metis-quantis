// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_finite_difference.hpp
#ifndef QUANTLIB_FINITE_DIFFERENCE_H
#define QUANTLIB_FINITE_DIFFERENCE_H

#include <ql/quantlib.hpp>
#include <vector>
#include <string>
#include <stdexcept>

namespace quant {

class QuantLibFiniteDifference {
public:
    enum class FDScheme {
        EXPLICIT_EULER,
        IMPLICIT_EULER,
        CRANK_NICOLSON,
        DOUGLAS
    };

    // Constructor
    QuantLibFiniteDifference(QuantLib::Size timeSteps = 100,
                             QuantLib::Size spaceSteps = 100,
                             QuantLib::Size dampingSteps = 0)
        : m_timeSteps(timeSteps),
          m_spaceSteps(spaceSteps),
          m_dampingSteps(dampingSteps) {}

    // European Option Pricing
    double priceEuropeanOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        FDScheme scheme = FDScheme::CRANK_NICOLSON);

    // American Option Pricing
    double priceAmericanOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        FDScheme scheme = FDScheme::CRANK_NICOLSON);

    // Barrier Option Pricing
    double priceBarrierOption(
        double spotPrice,
        double strikePrice,
        double barrierLevel,
        double rebate,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        QuantLib::Barrier::Type barrierType = QuantLib::Barrier::UpOut,
        FDScheme scheme = FDScheme::CRANK_NICOLSON);

    // Calculate option Greeks
    QuantLib::Greeks getGreeks(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        QuantLib::Exercise::Type exerciseType = QuantLib::Exercise::European,
        FDScheme scheme = FDScheme::CRANK_NICOLSON);

    // Get option values at all grid points for visualization
    std::vector<std::vector<double>> getOptionGrid(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        double xMin,
        double xMax,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        QuantLib::Exercise::Type exerciseType = QuantLib::Exercise::European,
        FDScheme scheme = FDScheme::CRANK_NICOLSON);

private:
    QuantLib::Size m_timeSteps;
    QuantLib::Size m_spaceSteps;
    QuantLib::Size m_dampingSteps;

    // Helper method to create finite difference pricing engine
    QuantLib::ext::shared_ptr<QuantLib::PricingEngine> createFDEngine(
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process,
        FDScheme scheme,
        QuantLib::Exercise::Type exerciseType = QuantLib::Exercise::European);
};

} // namespace quant

#endif // QUANTLIB_FINITE_DIFFERENCE_H