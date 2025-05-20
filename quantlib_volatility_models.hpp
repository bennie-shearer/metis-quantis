// Path: include/quantlib_volatility_models.hpp

#ifndef QUANTLIB_VOLATILITY_MODELS_H
#define QUANTLIB_VOLATILITY_MODELS_H

#include <ql/quantlib.hpp>
#include <vector>
#include <string>

namespace quant {

class QuantLibVolatilityModels {
public:
    // Constructor
    QuantLibVolatilityModels() {}

    // Create local volatility surface from market data
    QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> createLocalVolatilitySurface(
        const QuantLib::Date& referenceDate,
        const QuantLib::Calendar& calendar,
        const std::vector<QuantLib::Date>& dates,
        const std::vector<double>& strikes,
        const QuantLib::Matrix& volMatrix,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed());

    // Create a stochastic volatility surface from SABR parameters
    QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> createSABRVolatilitySurface(
        const QuantLib::Date& referenceDate,
        const QuantLib::Calendar& calendar,
        double forward,
        double alpha,
        double beta,
        double nu,
        double rho,
        const std::vector<double>& strikes,
        const std::vector<double>& maturities,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed());

    // Calculate implied volatility from Black-Scholes model
    double calculateBlackScholesImpliedVolatility(
        double optionPrice,
        double spotPrice,
        double strikePrice,
        double riskFreeRate,
        double dividendYield,
        double timeToMaturity,
        QuantLib::Option::Type optionType = QuantLib::Option::Call);

    // Build a volatility term structure from market data
    QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> buildVolatilityTermStructure(
        const QuantLib::Date& referenceDate,
        const std::vector<QuantLib::Period>& optionTenors,
        const std::vector<double>& strikes,
        const QuantLib::Matrix& vols,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed(),
        const QuantLib::Calendar& calendar = QuantLib::TARGET());

    // Calculate historical volatility from price series
    double calculateHistoricalVolatility(
        const std::vector<double>& prices,
        double annualizationFactor = 252.0,
        int lookbackWindow = 0);
};

} // namespace quant

#endif // QUANTLIB_VOLATILITY_MODELS_H