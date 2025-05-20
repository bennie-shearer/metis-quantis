// File: quantlib_exotic_options_api.cpp
#include "quantlib_exotic_options_api.hpp"
#include "quantlib_exotic_options.hpp"
#include "quantlib_monte_carlo.hpp"
#include "metis_json.hpp"

#include <ql/instruments/barrieroption.hpp>
#include <ql/instruments/asianoption.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/barrier/analyticbarrierengine.hpp>
#include <ql/pricingengines/barrier/mcbarrierengine.hpp>
#include <ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp>
#include <ql/pricingengines/asian/analytic_discr_geom_av_price.hpp>
#include <ql/pricingengines/asian/mc_discr_arith_av_price.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <cmath>
#include <iostream>

namespace quant {

simple_json::SimpleJSON handleExoticOptionsCalculation(const std::string& requestBody) {
    simple_json::SimpleJSON result = simple_json::SimpleJSON::object();

    try {
        // Parse the request JSON
        simple_json::SimpleJSON request = simple_json::SimpleJSON::fromString(requestBody);

        // Extract common parameters
        std::string optionType = request["optionType"].asString();
        double spotPrice = request["spotPrice"].asDouble();
        double strikePrice = request["strikePrice"].asDouble();
        double riskFreeRate = request["riskFreeRate"].asDouble() / 100.0; // Convert from percentage
        double dividendYield = request["dividendYield"].asDouble() / 100.0; // Convert from percentage
        double volatility = request["volatility"].asDouble() / 100.0; // Convert from percentage
        double timeToMaturity = request["timeToMaturity"].asDouble();
        std::string optionStyle = request["optionStyle"].asString();
        std::string pricingMethod = request["pricingMethod"].asString();

        // Create QuantLib option type
        QuantLib::Option::Type qlOptionType =
            (optionStyle == "call") ? QuantLib::Option::Call : QuantLib::Option::Put;

        // Set up basic market data
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::Date today = QuantLib::Date::todaysDate();
        QuantLib::Date expiryDate = today + static_cast<int>(timeToMaturity * 365);

        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> volatilityTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(today, QuantLib::TARGET(), volatility, dayCounter)));

        // Create Black-Scholes process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> bsProcess(
            new QuantLib::BlackScholesProcess(underlyingH, dividendTS,
                                            riskFreeTS, volatilityTS));

        // Create exotic option pricing object
        QuantLibExoticOptionsEnhanced exoticCalculator;

        // Price the option based on type
        double optionPrice = 0.0;
        double blackScholesPrice = 0.0;

        // Calculate equivalent vanilla option price
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(qlOptionType, strikePrice));
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(expiryDate));

        QuantLib::VanillaOption vanillaOption(payoff, exercise);
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> bsEngine(
            new QuantLib::AnalyticEuropeanEngine(bsProcess));
        vanillaOption.setPricingEngine(bsEngine);

        blackScholesPrice = vanillaOption.NPV();

        // Additional information section for specific option types
        simple_json::SimpleJSON additionalInfo = simple_json::SimpleJSON::object();

        // Price the exotic option based on option type
        if (optionType == "barrier") {
            // Extract barrier option parameters
            std::string barrierType = request["barrierType"].asString();
            double barrierLevel = request["barrierLevel"].asDouble();
            double rebate = request["rebate"].asDouble();

            // Convert string to QuantLib barrier type
            QuantLib::Barrier::Type qlBarrierType;
            if (barrierType == "up-out") {
                qlBarrierType = QuantLib::Barrier::UpOut;
            } else if (barrierType == "up-in") {
                qlBarrierType = QuantLib::Barrier::UpIn;
            } else if (barrierType == "down-out") {
                qlBarrierType = QuantLib::Barrier::DownOut;
            } else { // down-in
                qlBarrierType = QuantLib::Barrier::DownIn;
            }

            // Create barrier option and price it
            if (pricingMethod == "analytic") {
                optionPrice = exoticCalculator.priceBarrierOption(
                    spotPrice, strikePrice, barrierLevel, rebate,
                    riskFreeRate, dividendYield, volatility, timeToMaturity,
                    qlOptionType, qlBarrierType);
            } else if (pricingMethod == "monte-carlo") {
                optionPrice = exoticCalculator.priceBarrierOptionMonteCarlo(
                    spotPrice, strikePrice, barrierLevel, rebate,
                    riskFreeRate, dividendYield, volatility, timeToMaturity,
                    qlOptionType, qlBarrierType, 10000);
            } else { // finite-difference
                optionPrice = exoticCalculator.priceBarrierOptionFD(
                    spotPrice, strikePrice, barrierLevel, rebate,
                    riskFreeRate, dividendYield, volatility, timeToMaturity,
                    qlOptionType, qlBarrierType, 100, 100);
            }

            // Calculate additional info for barrier options
            // Hitting probability (simple approximation based on log-normal distribution)
            double hitProbability = 0.0;
            if (qlBarrierType == QuantLib::Barrier::UpOut || qlBarrierType == QuantLib::Barrier::UpIn) {
                double mu = (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * timeToMaturity;
                double sigma = volatility * std::sqrt(timeToMaturity);
                double d = (std::log(barrierLevel / spotPrice) - mu) / sigma;
                hitProbability = 0.5 * (1.0 + std::erf(d / std::sqrt(2.0)));
            } else { // Down barrier
                double mu = (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * timeToMaturity;
                double sigma = volatility * std::sqrt(timeToMaturity);
                double d = (std::log(spotPrice / barrierLevel) - mu) / sigma;
                hitProbability = 0.5 * (1.0 + std::erf(d / std::sqrt(2.0)));
            }

            additionalInfo["hitProbability"] = hitProbability;
            additionalInfo["expectedRebate"] = rebate * hitProbability * std::exp(-riskFreeRate * timeToMaturity);

        } else if (optionType == "asian") {
            // Extract Asian option parameters
            std::string averagingType = request["averagingType"].asString();
            int numObservations = request["numObservations"].asInteger();

            // Create averaging type
            QuantLib::Average::Type qlAveragingType =
                (averagingType == "arithmetic") ? QuantLib::Average::Arithmetic : QuantLib::Average::Geometric;

            // Price Asian option
            if (pricingMethod == "analytic" && qlAveragingType == QuantLib::Average::Geometric) {
                // Geometric average can be priced analytically
                optionPrice = exoticCalculator.priceAsianOption(
                    spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
                    qlOptionType, qlAveragingType, numObservations);
            } else {
                // Use Monte Carlo for arithmetic average
                optionPrice = exoticCalculator.priceAsianOptionMonteCarlo(
                    spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
                    qlOptionType, qlAveragingType, numObservations, 10000);
            }

            // Calculate additional info for Asian options
            // Expected average price
            double expectedAverage = spotPrice * std::exp((riskFreeRate - dividendYield) * timeToMaturity / 2.0);

            // Average price volatility (reduces with number of observations)
            double avgVolatility = volatility / std::sqrt(numObservations);

            additionalInfo["averagePrice"] = expectedAverage;
            additionalInfo["averageVolatility"] = avgVolatility;

        } else if (optionType == "lookback") {
            // Extract lookback option parameters
            std::string lookbackType = request["lookbackType"].asString();
            bool isFloatingStrike = (lookbackType == "floating-strike");

            // Price lookback option using Monte Carlo
            optionPrice = exoticCalculator.priceLookbackOption(
                spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
                qlOptionType, isFloatingStrike);

            // Calculate additional lookback info
            // Expected minimum/maximum
            double expectedMin = spotPrice * std::exp(-volatility * std::sqrt(timeToMaturity / M_PI));
            double expectedMax = spotPrice * std::exp(volatility * std::sqrt(timeToMaturity / M_PI));

            additionalInfo["expectedMinimum"] = expectedMin;
            additionalInfo["expectedMaximum"] = expectedMax;

        } else if (optionType == "digital") {
            // Extract digital option parameters
            double cashPayoff = request["cashPayoff"].asDouble();

            // Price digital option
            optionPrice = exoticCalculator.priceDigitalOption(
                spotPrice, strikePrice, volatility, riskFreeRate, dividendYield, timeToMaturity,
                qlOptionType, cashPayoff);

            // Calculate additional digital option info
            // Probability of exercise
            double d2 = (std::log(spotPrice / strikePrice) +
                       (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * timeToMaturity) /
                       (volatility * std::sqrt(timeToMaturity));

            double exerciseProbability = 0.0;
            if (qlOptionType == QuantLib::Option::Call) {
                exerciseProbability = 0.5 * (1.0 + std::erf(d2 / std::sqrt(2.0)));
            } else {
                exerciseProbability = 0.5 * (1.0 + std::erf(-d2 / std::sqrt(2.0)));
            }

            additionalInfo["exerciseProbability"] = exerciseProbability;

        } else if (optionType == "compound") {
            // Extract compound option parameters
            double motherExpiry = request["motherExpiry"].asDouble();
            std::string daughterStyle = request["daughterStyle"].asString();

            QuantLib::Option::Type qlDaughterType =
                (daughterStyle == "call") ? QuantLib::Option::Call : QuantLib::Option::Put;

            // Price compound option
            optionPrice = exoticCalculator.priceCompoundOption(
                spotPrice, strikePrice, strikePrice, volatility,
                riskFreeRate, dividendYield, motherExpiry, timeToMaturity,
                qlOptionType, qlDaughterType);

            // Calculate additional compound option info
            // Mother option value
            additionalInfo["motherOptionValue"] = optionPrice;

            // Daughter option expected value
            double daughterTimeToMaturity = timeToMaturity - motherExpiry;
            double expectedSpot = spotPrice * std::exp((riskFreeRate - dividendYield) * motherExpiry);

            // Create daughter payoff
            QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> daughterPayoff(
                new QuantLib::PlainVanillaPayoff(qlDaughterType, strikePrice));
            QuantLib::ext::shared_ptr<QuantLib::Exercise> daughterExercise(
                new QuantLib::EuropeanExercise(today + static_cast<int>(timeToMaturity * 365)));

            QuantLib::VanillaOption daughterOption(daughterPayoff, daughterExercise);

            // Create a temporary process for the daughter
            QuantLib::Handle<QuantLib::Quote> expectedSpotQuote(
                QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(expectedSpot)));

            QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> daughterProcess(
                new QuantLib::BlackScholesProcess(expectedSpotQuote, dividendTS,
                                                riskFreeTS, volatilityTS));

            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> daughterEngine(
                new QuantLib::AnalyticEuropeanEngine(daughterProcess));

            daughterOption.setPricingEngine(daughterEngine);
            additionalInfo["daughterOptionExpectedValue"] = daughterOption.NPV();
        }

        // Calculate option sensitivities (Greeks)
        QuantLib::Greeks greeks = exoticCalculator.calculateExoticGreeks(
            spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity,
            optionType, optionStyle, request);

        // Calculate intrinsic and time value
        double intrinsicValue = 0.0;
        if (qlOptionType == QuantLib::Option::Call) {
            intrinsicValue = std::max(spotPrice - strikePrice, 0.0);
        } else {
            intrinsicValue = std::max(strikePrice - spotPrice, 0.0);
        }

        double timeValue = optionPrice - intrinsicValue;

        // Generate payoff chart data
        simple_json::SimpleJSON payoffChart = simple_json::SimpleJSON::object();
        simple_json::SimpleJSON spotPrices = simple_json::SimpleJSON::array();
        simple_json::SimpleJSON payoffs = simple_json::SimpleJSON::array();
        simple_json::SimpleJSON currentValues = simple_json::SimpleJSON::array();
        simple_json::SimpleJSON vanillaPayoffs = simple_json::SimpleJSON::array();

        // Generate price range (50% to 150% of spot)
        double minPrice = spotPrice * 0.5;
        double maxPrice = spotPrice * 1.5;
        double step = (maxPrice - minPrice) / 50.0;

        for (double price = minPrice; price <= maxPrice; price += step) {
            // Store spot price
            spotPrices.pushBack(price);

            // Calculate payoff at expiry (simplified)
            double payoff = 0.0;
            if (optionType == "barrier") {
                std::string barrierType = request["barrierType"].asString();
                double barrierLevel = request["barrierLevel"].asDouble();
                double rebate = request["rebate"].asDouble();

                bool knockedOut = false;
                if (barrierType == "up-out" && price >= barrierLevel) {
                    knockedOut = true;
                } else if (barrierType == "down-out" && price <= barrierLevel) {
                    knockedOut = true;
                }

                if (knockedOut) {
                    payoff = rebate;
                } else {
                    if (qlOptionType == QuantLib::Option::Call) {
                        payoff = std::max(price - strikePrice, 0.0);
                    } else {
                        payoff = std::max(strikePrice - price, 0.0);
                    }
                }
            } else if (optionType == "digital") {
                double cashPayoff = request["cashPayoff"].asDouble();
                if (qlOptionType == QuantLib::Option::Call) {
                    payoff = (price > strikePrice) ? cashPayoff : 0.0;
                } else {
                    payoff = (price < strikePrice) ? cashPayoff : 0.0;
                }
            } else {
                // Simplified payoff - use vanilla option payoff
                if (qlOptionType == QuantLib::Option::Call) {
                    payoff = std::max(price - strikePrice, 0.0);
                } else {
                    payoff = std::max(strikePrice - price, 0.0);
                }
            }

            payoffs.pushBack(payoff);

            // Calculate vanilla option payoff
            double vanillaPayoff = 0.0;
            if (qlOptionType == QuantLib::Option::Call) {
                vanillaPayoff = std::max(price - strikePrice, 0.0);
            } else {
                vanillaPayoff = std::max(strikePrice - price, 0.0);
            }

            vanillaPayoffs.pushBack(vanillaPayoff);

            // Calculate current value (simplified)
            // Set up a temporary process with the new spot price
            QuantLib::Handle<QuantLib::Quote> tempSpot(
                QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(price)));

            QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> tempProcess(
                new QuantLib::BlackScholesProcess(tempSpot, dividendTS,
                                                riskFreeTS, volatilityTS));

            // Recalculate option price with new spot
            double currentValue = 0.0;

            // Simplified - use interpolation between intrinsic value and current price
            double distanceRatio = std::abs(price - spotPrice) / spotPrice;
            currentValue = optionPrice * (1.0 - distanceRatio) + payoff * distanceRatio;

            // Ensure current value is at least intrinsic value
            double tempIntrinsic = 0.0;
            if (qlOptionType == QuantLib::Option::Call) {
                tempIntrinsic = std::max(price - strikePrice, 0.0);
            } else {
                tempIntrinsic = std::max(strikePrice - price, 0.0);
            }

            currentValue = std::max(currentValue, tempIntrinsic);

            currentValues.pushBack(currentValue);
        }

        payoffChart["spotPrices"] = spotPrices;
        payoffChart["payoffs"] = payoffs;
        payoffChart["currentValues"] = currentValues;
        payoffChart["vanillaPayoffs"] = vanillaPayoffs;

        // Populate result
        result["optionPrice"] = optionPrice;
        result["blackScholesPrice"] = blackScholesPrice;
        result["intrinsicValue"] = intrinsicValue;
        result["timeValue"] = timeValue;

        // Add Greeks
        simple_json::SimpleJSON greeksJson = simple_json::SimpleJSON::object();
        greeksJson["delta"] = greeks.delta;
        greeksJson["gamma"] = greeks.gamma;
        greeksJson["vega"] = greeks.vega;
        greeksJson["theta"] = greeks.theta;
        greeksJson["rho"] = greeks.rho;

        result["greeks"] = greeksJson;
        result["additionalInfo"] = additionalInfo;
        result["payoffChart"] = payoffChart;

    } catch (std::exception& e) {
        result["error"] = std::string("Error calculating exotic option: ") + e.what();
        result["status"] = "error";
    }

    return result;
}

void registerExoticOptionsAPI(simple_web::SimpleWebServer& server) {
    server.post("/api/exotics/calculate", [](const std::string& request_content, simple_web::SimpleWebServer::ResponseHandler response_handler) {
        simple_json::SimpleJSON result = handleExoticOptionsCalculation(request_content);
        response_handler(200, result.dump());
    });
}

} // namespace quant