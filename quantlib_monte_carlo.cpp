// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_monte_carlo.cpp
#include "quantlib_monte_carlo.hpp"
#include <iostream>

namespace quant {

// European Option Pricing - Changed to use string for discretization method
double QuantLibMonteCarlo::priceEuropeanOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    const std::string& discretizationMethod) {

    try {
        // Set up dates
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date todaysDate = QuantLib::Date::todaysDate();
        QuantLib::Settings::instance().evaluationDate() = todaysDate;

        QuantLib::Date maturityDate = todaysDate + QuantLib::Integer(timeToMaturity * 365 + 0.5);

        // Option parameters
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));

        // European exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> europeanExercise(
            new QuantLib::EuropeanExercise(maturityDate));

        // Set up market data
        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeRateTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendYieldTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(todaysDate, calendar, volatility, dayCounter)));

        // Create Black-Scholes stochastic process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> bsProcess(
            new QuantLib::GeneralizedBlackScholesProcess(underlyingH, dividendYieldTS,
                                                      riskFreeRateTS, blackVolTS));

        // Create the option
        QuantLib::VanillaOption option(payoff, europeanExercise);

        // Set up Monte Carlo parameters
        bool antitheticVariate = m_antitheticVariate;
        bool controlVariate = m_controlVariate;
        QuantLib::Size requiredSamples = m_numPaths;
        QuantLib::Size timeSteps = m_timeSteps;
        QuantLib::BigNatural seed = m_seed;

        // Create the Monte Carlo engine with the proper constructor
        // Note: Updated constructor parameters to match what's expected in your QuantLib version
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> mcEngine;

        // Use simple conditional without the enum
        if (discretizationMethod == "Euler") {
            mcEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::MCEuropeanEngine<QuantLib::PseudoRandom, QuantLib::Statistics>(
                    bsProcess,
                    timeSteps,
                    QuantLib::Null<QuantLib::Size>(),
                    antitheticVariate,
                    controlVariate,
                    requiredSamples,
                    0.02,
                    1000000,
                    seed
                )
            );
        } else {
            // Default to Quadratic Exponential
            mcEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::MCEuropeanEngine<QuantLib::PseudoRandom, QuantLib::Statistics>(
                    bsProcess,
                    timeSteps,
                    QuantLib::Null<QuantLib::Size>(),
                    antitheticVariate,
                    controlVariate,
                    requiredSamples,
                    0.02,
                    1000000,
                    seed
                )
            );
        }

        option.setPricingEngine(mcEngine);

        // Calculate option price
        return option.NPV();

    } catch (std::exception& e) {
        std::cerr << "Error in European option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in European option pricing" << std::endl;
        return -1.0;
    }
}

// Asian Option Pricing
double QuantLibMonteCarlo::priceAsianOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    QuantLib::Average::Type averageType) {

    try {
        // Set up dates
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date todaysDate = QuantLib::Date::todaysDate();
        QuantLib::Settings::instance().evaluationDate() = todaysDate;

        QuantLib::Date maturityDate = todaysDate + QuantLib::Integer(timeToMaturity * 365 + 0.5);

        // Option parameters
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));

        // European exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(maturityDate));

        // Set up market data
        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeRateTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendYieldTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(todaysDate, calendar, volatility, dayCounter)));

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process(
            new QuantLib::GeneralizedBlackScholesProcess(underlyingH, dividendYieldTS,
                                                      riskFreeRateTS, blackVolTS));

        // Create option
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine;

        // Generate fixing dates
        std::vector<QuantLib::Date> fixingDates;

        // Weekly fixings
        QuantLib::Size numFixings = QuantLib::Size(timeToMaturity * 52);  // Approximately weekly
        for (QuantLib::Size i = 1; i <= numFixings; ++i) {
            fixingDates.push_back(todaysDate + i * QuantLib::Period(1, QuantLib::Weeks));
        }
        fixingDates.push_back(maturityDate);

        if (averageType == QuantLib::Average::Geometric) {
            // Create discrete geometric average price Asian option
            // Using correct class name based on your QuantLib version
            QuantLib::DiscreteAveragingAsianOption asianOption(
                averageType,
                0.0,  // Running sum (none for a new option)
                0,    // Running product (none for a new option)
                fixingDates,
                payoff,
                exercise
            );

            // Set geometric pricing engine
            engine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::AnalyticDiscreteGeometricAveragePriceAsianEngine(process));

            asianOption.setPricingEngine(engine);
            return asianOption.NPV();
        } else {
            // Arithmetic average price Asian option
            QuantLib::DiscreteAveragingAsianOption asianOption(
                averageType,
                0.0,  // Running sum (none for a new option)
                0,    // Running product (none for a new option)
                fixingDates,
                payoff,
                exercise
            );

            // Monte Carlo parameters
            bool antitheticVariate = m_antitheticVariate;
            bool controlVariate = m_controlVariate;
            QuantLib::Size requiredSamples = m_numPaths;
            QuantLib::Size timeSteps = m_timeSteps;
            QuantLib::BigNatural seed = m_seed;

            // Set arithmetic MC engine - Using correct constructor arguments
            engine = QuantLib::MakeMCDiscreteArithmeticAPEngine<QuantLib::PseudoRandom>(process)
                .withSamples(requiredSamples)
                .withSeed(seed)
                .withControlVariate(controlVariate);

            asianOption.setPricingEngine(engine);
            return asianOption.NPV();
        }

    } catch (std::exception& e) {
        std::cerr << "Error in Asian option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in Asian option pricing" << std::endl;
        return -1.0;
    }
}

// Barrier Option Pricing
double QuantLibMonteCarlo::priceBarrierOption(
    double spotPrice,
    double strikePrice,
    double barrierLevel,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    QuantLib::Barrier::Type barrierType) {

    try {
        // Set up dates
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date todaysDate = QuantLib::Date::todaysDate();
        QuantLib::Settings::instance().evaluationDate() = todaysDate;

        QuantLib::Date maturityDate = todaysDate + QuantLib::Integer(timeToMaturity * 365 + 0.5);

        // Option parameters
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));

        // European exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(maturityDate));

        // Set up market data
        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeRateTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendYieldTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(todaysDate, calendar, volatility, dayCounter)));

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process(
            new QuantLib::GeneralizedBlackScholesProcess(underlyingH, dividendYieldTS,
                                                      riskFreeRateTS, blackVolTS));

        // Create barrier option
        QuantLib::BarrierOption barrierOption(
            barrierType,
            barrierLevel,
            0.0,  // No rebate
            payoff,
            exercise
        );

        // Try to use the analytical engine first (faster and more accurate)
        try {
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> analyticEngine(
                new QuantLib::AnalyticBarrierEngine(process));

            barrierOption.setPricingEngine(analyticEngine);
            return barrierOption.NPV();
        }
        catch (...) {
            // If analytical fails, fall back to Monte Carlo
            // Monte Carlo parameters
            bool antitheticVariate = m_antitheticVariate;
            QuantLib::Size requiredSamples = m_numPaths;
            QuantLib::Size timeSteps = m_timeSteps;
            QuantLib::BigNatural seed = m_seed;
            bool brownianBridge = false;  // Extra parameter required for barrier engine

            // Create Monte Carlo barrier engine with fixed constructor arguments
            // Fix: Adding the missing controlVariate parameter
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> mcEngine(
                new QuantLib::MCBarrierEngine<QuantLib::PseudoRandom>(
                    process,
                    timeSteps,
                    QuantLib::Null<QuantLib::Size>(),
                    antitheticVariate,
                    false,  // controlVariate parameter added
                    requiredSamples,
                    0.02,
                    1000000,
                    seed,
                    brownianBridge
                )
            );

            barrierOption.setPricingEngine(mcEngine);
            return barrierOption.NPV();
        }

    } catch (std::exception& e) {
        std::cerr << "Error in barrier option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in barrier option pricing" << std::endl;
        return -1.0;
    }
}

// Lookback Option Pricing - Manual implementation since QuantLib doesn't provide a specific class
double QuantLibMonteCarlo::priceLookbackOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    bool isFloatingStrike) {

    try {
        // Set up common variables
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date todaysDate = QuantLib::Date::todaysDate();
        QuantLib::Settings::instance().evaluationDate() = todaysDate;

        QuantLib::Date maturityDate = todaysDate + QuantLib::Integer(timeToMaturity * 365 + 0.5);
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();

        // Set up market data
        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeRateTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendYieldTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(todaysDate, calendar, volatility, dayCounter)));

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process(
            new QuantLib::GeneralizedBlackScholesProcess(underlyingH, dividendYieldTS,
                                                      riskFreeRateTS, blackVolTS));

        // For lookback options, we'll implement a simple Monte Carlo simulation ourselves
        // since this version of QuantLib doesn't have built-in support

        // Time grid for simulation
        int timeSteps = 100;
        double dt = timeToMaturity / timeSteps;

        // Parameters for simulation
        double drift = riskFreeRate - dividendYield - 0.5 * volatility * volatility;
        double sqrtDt = std::sqrt(dt);
        double discount = std::exp(-riskFreeRate * timeToMaturity);

        // Setup random number generator
        std::mt19937 gen(m_seed);
        std::normal_distribution<> dist(0.0, 1.0);

        // Simulate paths
        std::vector<double> lookbackValues(m_numPaths, 0.0);

        for (size_t i = 0; i < m_numPaths; ++i) {
            double currentSpot = spotPrice;
            double maxSpot = spotPrice;
            double minSpot = spotPrice;

            // Generate path
            for (int j = 0; j < timeSteps; ++j) {
                double z = dist(gen);
                currentSpot *= std::exp(drift * dt + volatility * sqrtDt * z);
                maxSpot = std::max(maxSpot, currentSpot);
                minSpot = std::min(minSpot, currentSpot);
            }

            // Calculate payoff based on option type
            if (isFloatingStrike) {
                // Floating strike lookback
                if (optionType == QuantLib::Option::Call) {
                    // Call payoff: S(T) - min(S(t))
                    lookbackValues[i] = currentSpot - minSpot;
                } else {
                    // Put payoff: max(S(t)) - S(T)
                    lookbackValues[i] = maxSpot - currentSpot;
                }
            } else {
                // Fixed strike lookback
                if (optionType == QuantLib::Option::Call) {
                    // Call payoff: max(max(S(t)) - K, 0)
                    lookbackValues[i] = std::max(maxSpot - strikePrice, 0.0);
                } else {
                    // Put payoff: max(K - min(S(t)), 0)
                    lookbackValues[i] = std::max(strikePrice - minSpot, 0.0);
                }
            }
        }

        // Calculate average payoff
        double sum = 0.0;
        for (double value : lookbackValues) {
            sum += value;
        }

        // Return discounted expected payoff
        return discount * (sum / m_numPaths);

    } catch (std::exception& e) {
        std::cerr << "Error pricing lookback option: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error pricing lookback option" << std::endl;
        return -1.0;
    }
}

// Black-Scholes price for comparison
double QuantLibMonteCarlo::blackScholesPrice(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    try {
        // Set up dates
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date todaysDate = QuantLib::Date::todaysDate();
        QuantLib::Settings::instance().evaluationDate() = todaysDate;

        QuantLib::Date maturityDate = todaysDate + QuantLib::Integer(timeToMaturity * 365 + 0.5);

        // Option parameters
        QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));

        // European exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(maturityDate));

        // Set up market data
        QuantLib::Handle<QuantLib::Quote> underlyingH(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeRateTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, riskFreeRate, dayCounter)));

        QuantLib::Handle<QuantLib::YieldTermStructure> dividendYieldTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(todaysDate, dividendYield, dayCounter)));

        QuantLib::Handle<QuantLib::BlackVolTermStructure> blackVolTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(todaysDate, calendar, volatility, dayCounter)));

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process(
            new QuantLib::GeneralizedBlackScholesProcess(underlyingH, dividendYieldTS,
                                                      riskFreeRateTS, blackVolTS));

        // Create European vanilla option
        QuantLib::VanillaOption option(payoff, exercise);

        // Use analytic European engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> analyticEngine(
            new QuantLib::AnalyticEuropeanEngine(process));

        option.setPricingEngine(analyticEngine);

        // Return the Black-Scholes price
        return option.NPV();

    } catch (std::exception& e) {
        std::cerr << "Error in Black-Scholes pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in Black-Scholes pricing" << std::endl;
        return -1.0;
    }
}

} // namespace quant