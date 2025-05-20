// C:/Users/bshearer/CLionProjects/quant-boost/quantlib_commodity.cpp

#include "quantlib_commodity.hpp"
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <iostream>

namespace quant {

QuantLibCommodity::QuantLibCommodity() {
    // Constructor implementation
}

QuantLibCommodity::~QuantLibCommodity() {
    // Destructor implementation
}

double QuantLibCommodity::calculateCommodityForwardPrice(
    double spotPrice,
    double contractValue,
    const QuantLib::Date& evaluationDate,
    const QuantLib::Date& maturityDate,
    double riskFreeRate,
    double storageRate,
    double convenienceYield,
    const QuantLib::DayCounter& dayCounter) {

    try {
        // Set the evaluation date
        QuantLib::Settings::instance().evaluationDate() = evaluationDate;

        // Calculate time to maturity
        QuantLib::Time t = dayCounter.yearFraction(evaluationDate, maturityDate);

        // Apply the cost-of-carry model for commodity forward pricing
        double costOfCarry = riskFreeRate + storageRate - convenienceYield;
        double forwardPrice = spotPrice * std::exp(costOfCarry * t);

        // Calculate NPV
        double npv = (forwardPrice - contractValue) * std::exp(-riskFreeRate * t);

        return npv;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in calculateCommodityForwardPrice: " << e.what() << std::endl;
        return 0.0;
    }
}

double QuantLibCommodity::calculateCommodityOptionPrice(
    double spotPrice,
    double strikePrice,
    const QuantLib::Date& evaluationDate,
    const QuantLib::Date& maturityDate,
    double riskFreeRate,
    double volatility,
    const QuantLib::DayCounter& dayCounter,
    bool isCall) {

    try {
        // Set the evaluation date
        QuantLib::Settings::instance().evaluationDate() = evaluationDate;

        // Handle objects
        boost::shared_ptr<QuantLib::Quote> spot(new QuantLib::SimpleQuote(spotPrice));
        boost::shared_ptr<QuantLib::SimpleQuote> vol(new QuantLib::SimpleQuote(volatility));

        // Create volatility structure - fixed to use BlackConstantVol
        QuantLib::Calendar calendar = QuantLib::TARGET();
        boost::shared_ptr<QuantLib::BlackVolTermStructure> volTS(
            new QuantLib::BlackConstantVol(evaluationDate, calendar,
                                          QuantLib::Handle<QuantLib::Quote>(vol),
                                          dayCounter)
        );

        // Risk-free term structure - fixed to use FlatForward
        boost::shared_ptr<QuantLib::YieldTermStructure> rTS(
            new QuantLib::FlatForward(evaluationDate, riskFreeRate, dayCounter)
        );

        // Create the Black-Scholes process
        boost::shared_ptr<QuantLib::BlackScholesMertonProcess> bsProcess(
            new QuantLib::BlackScholesMertonProcess(
                QuantLib::Handle<QuantLib::Quote>(spot),
                QuantLib::Handle<QuantLib::YieldTermStructure>(rTS),
                QuantLib::Handle<QuantLib::YieldTermStructure>(rTS),  // Dividend yield = 0
                QuantLib::Handle<QuantLib::BlackVolTermStructure>(volTS)
            )
        );

        // Create the option type (call or put)
        boost::shared_ptr<QuantLib::StrikedTypePayoff> payoff(
            new QuantLib::PlainVanillaPayoff(
                isCall ? QuantLib::Option::Call : QuantLib::Option::Put,
                strikePrice
            )
        );

        // Create the exercise type (European)
        boost::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(maturityDate)
        );

        // Create the option
        QuantLib::VanillaOption option(payoff, exercise);

        // Set the pricing engine
        boost::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::AnalyticEuropeanEngine(bsProcess)
        );
        option.setPricingEngine(engine);

        // Calculate and return the NPV
        return option.NPV();
    }
    catch (const std::exception& e) {
        std::cerr << "Error in calculateCommodityOptionPrice: " << e.what() << std::endl;
        return 0.0;
    }
}

} // namespace quant