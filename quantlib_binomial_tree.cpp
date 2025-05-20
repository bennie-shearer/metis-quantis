// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_binomial_tree.cpp

#include <ql/quantlib.hpp>
#include "quantlib_binomial_tree.hpp"
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

namespace quant {

// Constructor
QuantLibBinomialTree::QuantLibBinomialTree(int steps)
    : m_steps(steps),
      m_dayCounter(QuantLib::Actual365Fixed()),
      m_calendar(QuantLib::TARGET()) {
}

// Helper method to create the appropriate binomial engine for use internally
QuantLib::ext::shared_ptr<QuantLib::PricingEngine> createCRREngine(
    QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> process,
    int steps) {

    return QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
        new QuantLib::BinomialVanillaEngine<QuantLib::CoxRossRubinstein>(process, steps));
}

// Implementation of the method defined in the header
double QuantLibBinomialTree::priceAmericanOption(
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
            new QuantLib::AmericanExercise(today, expiryDate));
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up the Black-Scholes process
        QuantLib::Handle<QuantLib::Quote> spot(
            QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(spotPrice)));
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, riskFreeRate, m_dayCounter)));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, dividendYield, m_dayCounter)));
        QuantLib::Handle<QuantLib::BlackVolTermStructure> volatilityTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(today, m_calendar, volatility, m_dayCounter)));

        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
            new QuantLib::BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, volatilityTS));

        // Create the binomial pricing engine (default CRR method)
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine =
            createCRREngine(bsProcess, m_steps);

        option.setPricingEngine(engine);

        // Calculate and return the option price
        return option.NPV();
    }
    catch (std::exception& e) {
        std::cerr << "Error in priceAmericanOption: " << e.what() << std::endl;
        return -1.0;
    }
    catch (...) {
        std::cerr << "Unknown error in priceAmericanOption" << std::endl;
        return -1.0;
    }
}

// Implementation of getGreeks method from the header
Greeks QuantLibBinomialTree::getGreeks(
    double spotPrice,
    double strikePrice,
    double riskFreeRate,
    double dividendYield,
    double volatility,
    double timeToMaturity,
    QuantLib::Option::Type optionType,
    QuantLib::Exercise::Type exerciseType) {

    Greeks greeks;

    try {
        // Set up dates
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date expiryDate = today + static_cast<int>(timeToMaturity * 365);

        // Set up the payoff
        QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(optionType, strikePrice));

        // Set up the exercise based on type
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise;

        if (exerciseType == QuantLib::Exercise::American) {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::AmericanExercise(today, expiryDate));
        } else {
            exercise = QuantLib::ext::shared_ptr<QuantLib::Exercise>(
                new QuantLib::EuropeanExercise(expiryDate));
        }

        // Create the option
        QuantLib::VanillaOption option(payoff, exercise);

        // Set up quotes for calculating finite differences
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> spotQuote(new QuantLib::SimpleQuote(spotPrice));
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> rateQuote(new QuantLib::SimpleQuote(riskFreeRate));
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> divQuote(new QuantLib::SimpleQuote(dividendYield));
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> volQuote(new QuantLib::SimpleQuote(volatility));

        // Create handles
        QuantLib::Handle<QuantLib::Quote> spot(spotQuote);
        QuantLib::Handle<QuantLib::YieldTermStructure> riskFreeTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, QuantLib::Handle<QuantLib::Quote>(rateQuote),
                                         m_dayCounter)));
        QuantLib::Handle<QuantLib::YieldTermStructure> dividendTS(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(today, QuantLib::Handle<QuantLib::Quote>(divQuote),
                                         m_dayCounter)));
        QuantLib::Handle<QuantLib::BlackVolTermStructure> volatilityTS(
            QuantLib::ext::shared_ptr<QuantLib::BlackVolTermStructure>(
                new QuantLib::BlackConstantVol(today, m_calendar,
                                             QuantLib::Handle<QuantLib::Quote>(volQuote),
                                             m_dayCounter)));

        // Create process
        QuantLib::ext::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
            new QuantLib::BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, volatilityTS));

        // Create pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine =
            createCRREngine(bsProcess, m_steps);

        option.setPricingEngine(engine);

        // Calculate Greeks using finite differences

        // Delta and Gamma
        double h = 0.01 * spotPrice;
        spotQuote->setValue(spotPrice + h);
        double pUp = option.NPV();

        spotQuote->setValue(spotPrice - h);
        double pDown = option.NPV();

        spotQuote->setValue(spotPrice);
        double p = option.NPV();

        greeks.delta = (pUp - pDown) / (2.0 * h);
        greeks.gamma = (pUp - 2.0 * p + pDown) / (h * h);

        // Theta (1-day)
        double t = 1.0 / 365.0;
        QuantLib::Settings::instance().evaluationDate() = today + 1;
        double pNext = option.NPV();
        QuantLib::Settings::instance().evaluationDate() = today;

        greeks.theta = (pNext - p) / t;

        // Vega (1% change)
        volQuote->setValue(volatility + 0.01);
        double pVolUp = option.NPV();

        greeks.vega = (pVolUp - p) / 0.01;

        // Rho (1% change)
        rateQuote->setValue(riskFreeRate + 0.01);
        double pRateUp = option.NPV();

        greeks.rho = (pRateUp - p) / 0.01;

        return greeks;
    }
    catch (std::exception& e) {
        std::cerr << "Error in getGreeks: " << e.what() << std::endl;
        greeks.delta = greeks.gamma = greeks.vega = greeks.theta = greeks.rho = 0.0;
        return greeks;
    }
    catch (...) {
        std::cerr << "Unknown error in getGreeks" << std::endl;
        greeks.delta = greeks.gamma = greeks.vega = greeks.theta = greeks.rho = 0.0;
        return greeks;
    }
}

} // namespace quant