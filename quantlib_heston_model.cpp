// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_heston_model.cpp
#include "quantlib_heston_model.hpp"
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/processes/hestonprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
// Use analytichestonengine.hpp instead of hestonintegral.hpp
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <iostream>

namespace quant {

double QuantLibHestonModel::priceEuropeanOption(
    double spot,
    double strike,
    double timeToMaturity,
    double riskFreeRate,
    const HestonParams& params,
    QuantLib::Option::Type optionType) {

    try {
        // Set up dates
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::Date today = QuantLib::Date::todaysDate();
        QuantLib::Date expiryDate = today + static_cast<QuantLib::Integer>(timeToMaturity * 365);

        // Set up quotes and curves
        QuantLib::Handle<QuantLib::Quote> spotQuote(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spot)));
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, riskFreeRate, dayCounter)));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, 0.0, dayCounter)));

        // Set up Heston process
        QuantLib::ext::shared_ptr<QuantLib::HestonProcess> hestonProcess(
            new QuantLib::HestonProcess(
                riskFreeTS,
                dividendTS,
                spotQuote,
                params.v0,    // Initial variance
                params.kappa, // Mean reversion speed
                params.theta, // Long-run variance
                params.sigma, // Volatility of variance
                params.rho    // Correlation
            ));

        // Create Heston model
        QuantLib::ext::shared_ptr<QuantLib::HestonModel> hestonModel(
            new QuantLib::HestonModel(hestonProcess));

        // Set up payoff
        QuantLib::ext::shared_ptr<QuantLib::PlainVanillaPayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strike));

        // Set up exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(expiryDate));

        // Create option
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up the pricing engine using AnalyticHestonEngine instead of HestonIntegralEngine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::AnalyticHestonEngine(hestonModel));
        option.setPricingEngine(engine);

        // Calculate option price
        return option.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error pricing European option with Heston model: " << e.what() << std::endl;
        return -1.0;
    }
}

std::vector<double> QuantLibHestonModel::priceEuropeanOptionGrid(
    double spot,
    const std::vector<double>& strikes,
    const std::vector<double>& maturities,
    double riskFreeRate,
    const HestonParams& params,
    QuantLib::Option::Type optionType) {

    try {
        // Set up dates
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::Date today = QuantLib::Date::todaysDate();

        // Set up quotes and curves
        QuantLib::Handle<QuantLib::Quote> spotQuote(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spot)));
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, riskFreeRate, dayCounter)));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, 0.0, dayCounter)));

        // Set up Heston process
        QuantLib::ext::shared_ptr<QuantLib::HestonProcess> hestonProcess(
            new QuantLib::HestonProcess(
                riskFreeTS,
                dividendTS,
                spotQuote,
                params.v0,    // Initial variance
                params.kappa, // Mean reversion speed
                params.theta, // Long-run variance
                params.sigma, // Volatility of variance
                params.rho    // Correlation
            ));

        // Create Heston model
        QuantLib::ext::shared_ptr<QuantLib::HestonModel> hestonModel(
            new QuantLib::HestonModel(hestonProcess));

        // Set up the pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::AnalyticHestonEngine(hestonModel));

        // Calculate prices for each strike and maturity
        std::vector<double> prices;
        prices.reserve(strikes.size() * maturities.size());

        for (double maturity : maturities) {
            QuantLib::Date expiryDate = today + static_cast<QuantLib::Integer>(maturity * 365);
            QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
                new QuantLib::EuropeanExercise(expiryDate));

            for (double strike : strikes) {
                QuantLib::ext::shared_ptr<QuantLib::PlainVanillaPayoff> payoff(
                    new QuantLib::PlainVanillaPayoff(optionType, strike));

                QuantLib::VanillaOption option(payoff, exercise);
                option.setPricingEngine(engine);

                prices.push_back(option.NPV());
            }
        }

        return prices;
    } catch (std::exception& e) {
        std::cerr << "Error pricing European option grid with Heston model: " << e.what() << std::endl;
        return std::vector<double>();
    }
}

std::vector<double> QuantLibHestonModel::calculateImpliedVolatilitySurface(
    double spot,
    const std::vector<double>& strikes,
    const std::vector<double>& maturities,
    double riskFreeRate,
    const HestonParams& params,
    QuantLib::Option::Type optionType) {

    try {
        // First calculate option prices
        std::vector<double> prices = priceEuropeanOptionGrid(
            spot, strikes, maturities, riskFreeRate, params, optionType);

        if (prices.empty()) {
            return std::vector<double>();
        }

        // Set up dates
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::Date today = QuantLib::Date::todaysDate();

        // Set up quotes and curves for Black-Scholes model
        QuantLib::Handle<QuantLib::Quote> spotQuote(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spot)));
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, riskFreeRate, dayCounter)));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, 0.0, dayCounter)));

        // Calculate implied volatilities
        std::vector<double> impliedVols;
        impliedVols.reserve(prices.size());

        size_t priceIndex = 0;
        for (double maturity : maturities) {
            QuantLib::Date expiryDate = today + static_cast<QuantLib::Integer>(maturity * 365);
            QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
                new QuantLib::EuropeanExercise(expiryDate));

            for (double strike : strikes) {
                // Get the Heston price
                double hestonPrice = prices[priceIndex++];

                // Create Black-Scholes process with dummy volatility (0.2)
                QuantLib::ext::shared_ptr<QuantLib::Quote> volQuote(new QuantLib::SimpleQuote(0.2));
                QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
                    QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                        new QuantLib::BlackConstantVol(today, QuantLib::NullCalendar(),
                                                     QuantLib::Handle<QuantLib::Quote>(volQuote),
                                                     dayCounter)));

                QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
                    new QuantLib::BlackScholesMertonProcess(
                        spotQuote, dividendTS, riskFreeTS, blackVolTS));

                // Create vanilla option
                QuantLib::ext::shared_ptr<QuantLib::PlainVanillaPayoff> payoff(
                    new QuantLib::PlainVanillaPayoff(optionType, strike));
                QuantLib::VanillaOption option(payoff, exercise);

                // Set up Black-Scholes pricing engine
                option.setPricingEngine(QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::AnalyticEuropeanEngine(bsProcess)));

                // Calculate implied volatility
                try {
                    double impliedVol = option.impliedVolatility(hestonPrice, bsProcess);
                    impliedVols.push_back(impliedVol);
                } catch (std::exception& e) {
                    // If implied vol calculation fails, use a default value
                    impliedVols.push_back(-1.0);
                }
            }
        }

        return impliedVols;
    } catch (std::exception& e) {
        std::cerr << "Error calculating implied volatility surface: " << e.what() << std::endl;
        return std::vector<double>();
    }
}

double QuantLibHestonModel::calculateHestonCalibrationError(
    const std::vector<double>& marketStrikes,
    const std::vector<double>& marketMaturities,
    const std::vector<double>& marketVols,
    double spot,
    double riskFreeRate,
    const HestonParams& params) {

    try {
        // Calculate model implied volatilities
        std::vector<double> modelVols = calculateImpliedVolatilitySurface(
            spot, marketStrikes, marketMaturities, riskFreeRate, params, QuantLib::Option::Call);

        if (modelVols.size() != marketVols.size()) {
            return std::numeric_limits<double>::max();
        }

        // Calculate mean squared error
        double totalError = 0.0;
        for (size_t i = 0; i < modelVols.size(); ++i) {
            if (modelVols[i] > 0) { // Only use valid implied vols
                double error = modelVols[i] - marketVols[i];
                totalError += error * error;
            }
        }

        return std::sqrt(totalError / modelVols.size());
    } catch (std::exception& e) {
        std::cerr << "Error calculating Heston calibration error: " << e.what() << std::endl;
        return std::numeric_limits<double>::max();
    }
}

} // namespace quant