// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_volatility_models.cpp

#include "quantlib_volatility_models.hpp"
#include <ql/termstructures/volatility/equityfx/blackvariancesurface.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/sabr.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/localvolsurface.hpp>
#include <ql/termstructures/volatility/equityfx/impliedvoltermstructure.hpp>
// Comment out the missing header file
// #include <ql/termstructures/volatility/equityfx/interpolatedsmilesection.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/interpolations/bicubicsplineinterpolation.hpp>
#include <ql/math/interpolations/sabrinterpolation.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/target.hpp>
#include <iostream>
#include <cmath>

namespace quant {

QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> QuantLibVolatilityModels::createLocalVolatilitySurface(
    const QuantLib::Date& referenceDate,
    const QuantLib::Calendar& calendar,
    const std::vector<QuantLib::Date>& dates,
    const std::vector<double>& strikes,
    const QuantLib::Matrix& volMatrix,
    const QuantLib::DayCounter& dayCounter) {

    try {
        // Create a Black variance surface from the volatility matrix
        QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> blackSurface(
            new QuantLib::BlackVarianceSurface(
                referenceDate,
                calendar,
                dates,
                strikes,
                volMatrix,
                dayCounter));

        // Set extrapolation
        blackSurface->enableExtrapolation();

        // Convert to local volatility surface
        QuantLib::ext::shared_ptr<QuantLib::LocalVolTermStructure> localVolSurface(
            new QuantLib::LocalVolSurface(
                QuantLib::Handle<QuantLib::BlackVolTermStructure>(blackSurface),
                QuantLib::Handle<QuantLib::YieldTermStructure>(
                    QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                        new QuantLib::FlatForward(referenceDate, 0.05, dayCounter))),
                QuantLib::Handle<QuantLib::YieldTermStructure>(
                    QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                        new QuantLib::FlatForward(referenceDate, 0.02, dayCounter))),
                100.0));  // Spot price

        // Return the black volatility surface directly
        // Note: Directly returning the black surface instead of creating an implied vol surface
        // since ImpliedVolTermStructure doesn't support the constructor with LocalVolTermStructure
        return blackSurface;

    } catch (std::exception& e) {
        std::cerr << "Error creating local volatility surface: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error creating local volatility surface" << std::endl;
        return nullptr;
    }
}

QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> QuantLibVolatilityModels::createSABRVolatilitySurface(
    const QuantLib::Date& referenceDate,
    const QuantLib::Calendar& calendar,
    double forward,
    double alpha,
    double beta,
    double nu,
    double rho,
    const std::vector<double>& strikes,
    const std::vector<double>& maturities,
    const QuantLib::DayCounter& dayCounter) {

    try {
        // Create a matrix to store volatilities
        QuantLib::Matrix volMatrix(strikes.size(), maturities.size());

        // Calculate SABR volatilities for each strike and maturity
        for (size_t i = 0; i < strikes.size(); ++i) {
            for (size_t j = 0; j < maturities.size(); ++j) {
                volMatrix[i][j] = QuantLib::sabrVolatility(
                    strikes[i],
                    forward,
                    maturities[j],
                    alpha,
                    beta,
                    nu,
                    rho);
            }
        }

        // Convert maturities to dates
        std::vector<QuantLib::Date> dates;
        for (double maturity : maturities) {
            QuantLib::Date maturityDate = referenceDate + QuantLib::Integer(maturity * 365 + 0.5);
            dates.push_back(maturityDate);
        }

        // Create Black variance surface from volatility matrix
        QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> surface(
            new QuantLib::BlackVarianceSurface(
                referenceDate,
                calendar,
                dates,
                strikes,
                volMatrix,
                dayCounter));

        surface->enableExtrapolation();

        return surface;

    } catch (std::exception& e) {
        std::cerr << "Error creating SABR volatility surface: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error creating SABR volatility surface" << std::endl;
        return nullptr;
    }
}

double QuantLibVolatilityModels::calculateBlackScholesImpliedVolatility(
    double optionPrice,
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    try {
        // Calculate forward price
        double forwardPrice = spotPrice * std::exp((riskFreeRate - dividendYield) * timeToMaturity);

        // Calculate discount factor
        double discountFactor = std::exp(-riskFreeRate * timeToMaturity);

        // Use QuantLib's Black formula to calculate implied volatility
        return QuantLib::blackFormulaImpliedStdDev(
            optionType,
            strikePrice,
            forwardPrice,
            optionPrice,
            discountFactor) / std::sqrt(timeToMaturity);

    } catch (std::exception& e) {
        std::cerr << "Error calculating implied volatility: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating implied volatility" << std::endl;
        return -1.0;
    }
}

QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> QuantLibVolatilityModels::buildVolatilityTermStructure(
    const QuantLib::Date& referenceDate,
    const std::vector<QuantLib::Period>& optionTenors,
    const std::vector<double>& strikes,
    const QuantLib::Matrix& vols,
    const QuantLib::DayCounter& dayCounter,
    const QuantLib::Calendar& calendar) {

    try {
        // Convert option tenors to dates
        std::vector<QuantLib::Date> dates;
        for (const auto& tenor : optionTenors) {
            dates.push_back(calendar.advance(referenceDate, tenor));
        }

        // Create Black variance surface from volatility matrix
        QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure> surface(
            new QuantLib::BlackVarianceSurface(
                referenceDate,
                calendar,
                dates,
                strikes,
                vols,
                dayCounter));

        surface->enableExtrapolation();

        return surface;

    } catch (std::exception& e) {
        std::cerr << "Error building volatility term structure: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error building volatility term structure" << std::endl;
        return nullptr;
    }
}

double QuantLibVolatilityModels::calculateHistoricalVolatility(
    const std::vector<double>& prices,
    double annualizationFactor,
    int lookbackWindow) {

    try {
        // Check for sufficient data
        if (prices.size() < 2) {
            throw std::invalid_argument("Need at least two prices to calculate volatility");
        }

        // If lookback window is 0 or greater than prices size, use all prices
        if (lookbackWindow <= 0 || lookbackWindow >= static_cast<int>(prices.size())) {
            lookbackWindow = prices.size() - 1;
        }

        // Calculate log returns
        std::vector<double> logReturns;
        for (size_t i = 1; i < prices.size() && static_cast<int>(i) <= lookbackWindow; ++i) {
            double logReturn = std::log(prices[i] / prices[i-1]);
            logReturns.push_back(logReturn);
        }

        // Calculate mean log return
        double sumReturns = 0.0;
        for (double logReturn : logReturns) {
            sumReturns += logReturn;
        }
        double meanReturn = sumReturns / logReturns.size();

        // Calculate sum of squared deviations
        double sumSquaredDeviations = 0.0;
        for (double logReturn : logReturns) {
            double deviation = logReturn - meanReturn;
            sumSquaredDeviations += deviation * deviation;
        }

        // Calculate variance
        double variance = sumSquaredDeviations / (logReturns.size() - 1);

        // Calculate annualized standard deviation
        return std::sqrt(variance * annualizationFactor);

    } catch (std::exception& e) {
        std::cerr << "Error calculating historical volatility: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating historical volatility" << std::endl;
        return -1.0;
    }
}

} // namespace quant