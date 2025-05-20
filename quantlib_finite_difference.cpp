// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_finite_difference.cpp
#include "quantlib_finite_difference.hpp"
#include <iostream>

namespace quant {

// Helper method to create the finite difference pricing engine
QuantLib::ext::shared_ptr<QuantLib::PricingEngine> QuantLibFiniteDifference::createFDEngine(
    QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process,
    FDScheme scheme,
    QuantLib::Exercise::Type exerciseType) {

    // For European options
    if (exerciseType == QuantLib::Exercise::European) {
        switch (scheme) {
            case FDScheme::EXPLICIT_EULER:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid (formerly called dampingSteps)
                        QuantLib::FdmSchemeDesc::ExplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::IMPLICIT_EULER:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid
                        QuantLib::FdmSchemeDesc::ImplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::CRANK_NICOLSON:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid
                        QuantLib::FdmSchemeDesc::CrankNicolson(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::DOUGLAS:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            default:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
        }
    }
    else if (exerciseType == QuantLib::Exercise::American) {
        // For American options
        switch (scheme) {
            case FDScheme::EXPLICIT_EULER:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid
                        QuantLib::FdmSchemeDesc::ExplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::IMPLICIT_EULER:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid
                        QuantLib::FdmSchemeDesc::ImplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::CRANK_NICOLSON:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps, // tGrid
                        QuantLib::FdmSchemeDesc::CrankNicolson(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            case FDScheme::DOUGLAS:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
            default:
                return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesVanillaEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
        }
    }

    // Default to European with Crank-Nicolson
    return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
        new QuantLib::FdBlackScholesVanillaEngine(
            process,
            m_timeSteps,
            m_spaceSteps,
            m_dampingSteps,
            QuantLib::FdmSchemeDesc::CrankNicolson(),
            true,  // localVol
            0.0    // illegalLocalVolOverwrite
        )
    );
}

// European Option Pricing - implementation remains the same
double QuantLibFiniteDifference::priceEuropeanOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    FDScheme scheme) {
    // Existing implementation is fine
    // This method will use the fixed createFDEngine

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

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process(
            new QuantLib::BlackScholesMertonProcess(underlyingH, dividendYieldTS,
                                                 riskFreeRateTS, blackVolTS));

        // Create vanilla option
        QuantLib::VanillaOption option(payoff, europeanExercise);

        // Create and set the finite difference pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> fdEngine = createFDEngine(
            process, scheme, QuantLib::Exercise::European);
        option.setPricingEngine(fdEngine);

        // Calculate result
        return option.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error in European option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in European option pricing" << std::endl;
        return -1.0;
    }
}

// American Option Pricing
double QuantLibFiniteDifference::priceAmericanOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    FDScheme scheme) {

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

        // American exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> americanExercise(
            new QuantLib::AmericanExercise(todaysDate, maturityDate));

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
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process(
            new QuantLib::BlackScholesMertonProcess(underlyingH, dividendYieldTS,
                                                 riskFreeRateTS, blackVolTS));

        // Create vanilla option
        QuantLib::VanillaOption option(payoff, americanExercise);

        // Create and set the finite difference pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> fdEngine = createFDEngine(
            process, scheme, QuantLib::Exercise::American);
        option.setPricingEngine(fdEngine);

        // Calculate result
        return option.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error in American option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in American option pricing" << std::endl;
        return -1.0;
    }
}

// Barrier Option Pricing
double QuantLibFiniteDifference::priceBarrierOption(
    double spotPrice,
    double strikePrice,
    double barrierLevel,
    double rebate,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    QuantLib::Barrier::Type barrierType,
    FDScheme scheme) {

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

        // European exercise for barrier option
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

        // Create stochastic process
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process(
            new QuantLib::BlackScholesMertonProcess(underlyingH, dividendYieldTS,
                                                 riskFreeRateTS, blackVolTS));

        // Create barrier option
        QuantLib::BarrierOption barrierOption(barrierType, barrierLevel, rebate,
                                          payoff, europeanExercise);

        // Create the FD barrier engine with the correct signature
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> fdEngine;

        switch (scheme) {
            case FDScheme::EXPLICIT_EULER:
                fdEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesBarrierEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::ExplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
                break;
            case FDScheme::IMPLICIT_EULER:
                fdEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesBarrierEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::ImplicitEuler(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
                break;
            case FDScheme::CRANK_NICOLSON:
                fdEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesBarrierEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::CrankNicolson(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
                break;
            case FDScheme::DOUGLAS:
                fdEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesBarrierEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
                break;
            default:
                fdEngine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                    new QuantLib::FdBlackScholesBarrierEngine(
                        process,
                        m_timeSteps,
                        m_spaceSteps,
                        m_dampingSteps,
                        QuantLib::FdmSchemeDesc::Douglas(),
                        true,  // localVol
                        0.0    // illegalLocalVolOverwrite
                    )
                );
        }

        barrierOption.setPricingEngine(fdEngine);

        // Calculate result
        return barrierOption.NPV();

    } catch (std::exception& e) {
        std::cerr << "Error in barrier option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in barrier option pricing" << std::endl;
        return -1.0;
    }
}

// Calculate option Greeks
QuantLib::Greeks QuantLibFiniteDifference::getGreeks(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    QuantLib::Exercise::Type exerciseType,
    FDScheme scheme) {

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

        // Create exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise;
        if (exerciseType == QuantLib::Exercise::European) {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::EuropeanExercise(maturityDate));
        } else {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::AmericanExercise(todaysDate, maturityDate));
        }

        // Create vanilla option
        QuantLib::VanillaOption option(payoff, exercise);

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
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process(
            new QuantLib::BlackScholesMertonProcess(underlyingH, dividendYieldTS,
                                                 riskFreeRateTS, blackVolTS));

        // Create and set the finite difference pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> fdEngine = createFDEngine(
            process, scheme, exerciseType);
        option.setPricingEngine(fdEngine);

        // Calculate Greeks
        QuantLib::Greeks greeks;
        greeks.delta = option.delta();
        greeks.gamma = option.gamma();
        greeks.theta = option.theta();
        greeks.vega = option.vega();
        greeks.rho = option.rho();

        return greeks;
    } catch (std::exception& e) {
        std::cerr << "Error calculating Greeks: " << e.what() << std::endl;
        return QuantLib::Greeks();
    } catch (...) {
        std::cerr << "Unknown error calculating Greeks" << std::endl;
        return QuantLib::Greeks();
    }
}

// Get option values at all grid points for visualization
std::vector<std::vector<double>> QuantLibFiniteDifference::getOptionGrid(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    double xMin,
    double xMax,
    QuantLib::Option::Type optionType,
    QuantLib::Exercise::Type exerciseType,
    FDScheme scheme) {

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

        // Create exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise;
        if (exerciseType == QuantLib::Exercise::European) {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::EuropeanExercise(maturityDate));
        } else {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::AmericanExercise(todaysDate, maturityDate));
        }

        // Create vanilla option
        QuantLib::VanillaOption option(payoff, exercise);

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
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process(
            new QuantLib::BlackScholesMertonProcess(underlyingH, dividendYieldTS,
                                                 riskFreeRateTS, blackVolTS));

        // Create and set the finite difference pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> fdEngine = createFDEngine(
            process, scheme, exerciseType);
        option.setPricingEngine(fdEngine);

        // Create grid of option values
        std::vector<std::vector<double>> optionGrid(m_timeSteps + 1);

        // Sample times
        std::vector<double> times;
        double dt = timeToMaturity / m_timeSteps;
        for (QuantLib::Size i = 0; i <= m_timeSteps; ++i) {
            times.push_back(i * dt);
        }

        // Sample stock prices
        std::vector<double> stockPrices;
        double dx = (xMax - xMin) / m_spaceSteps;
        for (QuantLib::Size j = 0; j <= m_spaceSteps; ++j) {
            stockPrices.push_back(xMin + j * dx);
        }

        // For each time point
        for (QuantLib::Size i = 0; i <= m_timeSteps; ++i) {
            optionGrid[i].resize(m_spaceSteps + 1);

            // For each stock price
            for (QuantLib::Size j = 0; j <= m_spaceSteps; ++j) {
                double S = stockPrices[j];

                // Set the current spot price
                QuantLib::ext::dynamic_pointer_cast<QuantLib::SimpleQuote>(
                    underlyingH.currentLink())->setValue(S);

                // Create a new vanilla option at the current spot price
                QuantLib::VanillaOption currentOption(payoff, exercise);
                currentOption.setPricingEngine(fdEngine);

                // Calculate the option value
                optionGrid[i][j] = currentOption.NPV();
            }
        }

        return optionGrid;
    } catch (std::exception& e) {
        std::cerr << "Error calculating option grid: " << e.what() << std::endl;
        return std::vector<std::vector<double>>();
    } catch (...) {
        std::cerr << "Unknown error calculating option grid" << std::endl;
        return std::vector<std::vector<double>>();
    }
}

} // namespace quant