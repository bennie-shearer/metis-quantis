// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_option_calculator.hpp

#ifndef QUANT_BOOST_QUANTLIB_OPTION_CALCULATOR_HPP
#define QUANT_BOOST_QUANTLIB_OPTION_CALCULATOR_HPP

#include <ql/quantlib.hpp>
#include <string>
#include <vector>
#include <map>

namespace quant {

// Define a Greeks structure to hold option sensitivities
struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

class QuantLibOptionCalculator {
public:
    QuantLibOptionCalculator();
    ~QuantLibOptionCalculator() = default;

    // Basic Black-Scholes pricing
    double calculateBlackScholesPrice(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    // European option pricing
    double calculateEuropeanOption(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    // Calculate all Greeks at once
    Greeks calculateGreeks(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    // Individual Greeks calculations
    double calculateDelta(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    double calculateGamma(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    double calculateVega(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    double calculateTheta(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    double calculateRho(
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double volatility,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

    // Implied volatility calculation
    double calculateImpliedVolatility(
        double optionPrice,
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double timeToMaturity,
        QuantLib::Option::Type optionType);

private:
    QuantLib::DayCounter m_dayCounter;
    QuantLib::Calendar m_calendar;
};

} // namespace quant

#endif // QUANT_BOOST_QUANTLIB_OPTION_CALCULATOR_HPP