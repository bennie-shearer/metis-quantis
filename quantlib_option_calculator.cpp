// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_option_calculator.cpp

#include "quantlib_option_calculator.hpp"
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/target.hpp>

namespace quant {

QuantLibOptionCalculator::QuantLibOptionCalculator()
    : m_dayCounter(QuantLib::Actual365Fixed()),
      m_calendar(QuantLib::TARGET()) {
}

double QuantLibOptionCalculator::calculateBlackScholesPrice(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    try {
        // Set up dates
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date expiryDate = today + static_cast<int>(timeToMaturity * 365);

        // Set up the option
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(expiryDate));
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up the Black-Scholes process
        QuantLib::Handle<QuantLib::Quote> spot(QuantLib::ext::make_shared<QuantLib::SimpleQuote>(spotPrice));
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, riskFreeRate, m_dayCounter));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, dividendYield, m_dayCounter));
        QuantLib::Handle<QuantLib::BlackVolTermStructure> volatilityTS(
            QuantLib::ext::make_shared<QuantLib::BlackConstantVol>(today, m_calendar, volatility, m_dayCounter));

        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
            new QuantLib::BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, volatilityTS));

        // Set up the pricing engine
        option.setPricingEngine(QuantLib::ext::make_shared<QuantLib::AnalyticEuropeanEngine>(bsProcess));

        // Calculate and return the option price
        return option.NPV();
    }
    catch (std::exception& e) {
        std::cerr << "Error in calculateBlackScholesPrice: " << e.what() << std::endl;
        return -1.0;
    }
    catch (...) {
        std::cerr << "Unknown error in calculateBlackScholesPrice" << std::endl;
        return -1.0;
    }
}

double QuantLibOptionCalculator::calculateEuropeanOption(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    // Use the Black-Scholes calculator for European options
    return calculateBlackScholesPrice(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType);
}

Greeks QuantLibOptionCalculator::calculateGreeks(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    Greeks greeks;

    try {
        // Set up dates
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date expiryDate = today + static_cast<int>(timeToMaturity * 365);

        // Set up the option
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(expiryDate));
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up the Black-Scholes process
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> spotQuote(new QuantLib::SimpleQuote(spotPrice));
        QuantLib::Handle<QuantLib::Quote> spot(spotQuote);

        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, riskFreeRate, m_dayCounter));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, dividendYield, m_dayCounter));

        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> volQuote(new QuantLib::SimpleQuote(volatility));
        QuantLib::Handle<QuantLib::BlackVolTermStructure> volatilityTS(
            QuantLib::ext::make_shared<QuantLib::BlackConstantVol>(today, m_calendar,
                QuantLib::Handle<QuantLib::Quote>(volQuote), m_dayCounter));

        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
            new QuantLib::BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, volatilityTS));

        // Set up the pricing engine
        option.setPricingEngine(QuantLib::ext::make_shared<QuantLib::AnalyticEuropeanEngine>(bsProcess));

        // Calculate greeks
        greeks.delta = option.delta();
        greeks.gamma = option.gamma();
        greeks.vega = option.vega() / 100.0; // QuantLib returns vega per percentage point
        greeks.theta = option.theta() / 365.0; // QuantLib returns theta per year
        greeks.rho = option.rho() / 100.0; // QuantLib returns rho per percentage point

        return greeks;
    }
    catch (std::exception& e) {
        std::cerr << "Error in calculateGreeks: " << e.what() << std::endl;
        greeks.delta = greeks.gamma = greeks.vega = greeks.theta = greeks.rho = 0.0;
        return greeks;
    }
    catch (...) {
        std::cerr << "Unknown error in calculateGreeks" << std::endl;
        greeks.delta = greeks.gamma = greeks.vega = greeks.theta = greeks.rho = 0.0;
        return greeks;
    }
}

// Implement individual Greeks calculations
double QuantLibOptionCalculator::calculateDelta(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    return calculateGreeks(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType).delta;
}

double QuantLibOptionCalculator::calculateGamma(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    return calculateGreeks(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType).gamma;
}

double QuantLibOptionCalculator::calculateVega(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    return calculateGreeks(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType).vega;
}

double QuantLibOptionCalculator::calculateTheta(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    return calculateGreeks(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType).theta;
}

double QuantLibOptionCalculator::calculateRho(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    return calculateGreeks(
        spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity, optionType).rho;
}

double QuantLibOptionCalculator::calculateImpliedVolatility(
    double optionPrice,
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double timeToMaturity,
    QuantLib::Option::Type optionType) {

    try {
        // Set up dates
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date expiryDate = today + static_cast<int>(timeToMaturity * 365);

        // Set up the option
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(expiryDate));
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up yield term structures
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, riskFreeRate, m_dayCounter));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::make_shared<QuantLib::FlatForward>(today, dividendYield, m_dayCounter));

        // Calculate implied volatility
        QuantLib::ext::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process(
            new QuantLib::BlackScholesMertonProcess(
                QuantLib::Handle<QuantLib::Quote>(
                    QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice))),
                dividendTS,
                riskFreeTS,
                QuantLib::Handle<QuantLib::BlackVolTermStructure>(
                    QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                        new QuantLib::BlackConstantVol(today, m_calendar, 0.20, m_dayCounter)
                    )
                )
            )
        );

        option.setPricingEngine(
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::AnalyticEuropeanEngine(process)
            )
        );

        return option.impliedVolatility(
            optionPrice,
            process,
            1.0e-6, // accuracy
            1000,   // max evaluations
            1.0e-7, // min vol
            4.0     // max vol
        );
    }
    catch (std::exception& e) {
        std::cerr << "Error in calculateImpliedVolatility: " << e.what() << std::endl;
        return -1.0;
    }
    catch (...) {
        std::cerr << "Unknown error in calculateImpliedVolatility" << std::endl;
        return -1.0;
    }
}

} // namespace quant