// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_monte_carlo.hpp
#ifndef QUANTLIB_MONTE_CARLO_H
#define QUANTLIB_MONTE_CARLO_H

#include <ql/quantlib.hpp>

namespace quant {

class QuantLibMonteCarlo {
public:
    // Constructor
    QuantLibMonteCarlo(QuantLib::Size numPaths = 10000,
                      QuantLib::Size timeSteps = 252,
                      bool antitheticVariate = true,
                      bool controlVariate = false,
                      bool brownianBridge = false,
                      QuantLib::BigNatural seed = 42)
        : m_numPaths(numPaths),
          m_timeSteps(timeSteps),
          m_antitheticVariate(antitheticVariate),
          m_controlVariate(controlVariate),
          m_brownianBridge(brownianBridge),
          m_seed(seed) {}

    // European Option Pricing - Updated to use a string for discretization type
    double priceEuropeanOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        const std::string& discretizationMethod = "Euler");

    // Asian Option Pricing
    double priceAsianOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        QuantLib::Average::Type averageType = QuantLib::Average::Arithmetic);

    // Barrier Option Pricing
    double priceBarrierOption(
        double spotPrice,
        double strikePrice,
        double barrierLevel,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        QuantLib::Barrier::Type barrierType = QuantLib::Barrier::UpOut);

    // Lookback Option Pricing
    double priceLookbackOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call,
        bool isFloatingStrike = false);

    // Black-Scholes price for comparison
    double blackScholesPrice(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call);

private:
    QuantLib::Size m_numPaths;
    QuantLib::Size m_timeSteps;
    bool m_antitheticVariate;
    bool m_controlVariate;
    bool m_brownianBridge;
    QuantLib::BigNatural m_seed;
};

} // namespace quant

#endif // QUANTLIB_MONTE_CARLO_H