// C:/Users/bshearer/CLionProjects/quant-boost/quantlib_bonds.cpp

#include "quantlib_bonds.hpp"
#include <ql/instruments/bonds/amortizingfixedratebond.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <iostream>

namespace quant {

QuantLibBondsEnhanced::QuantLibBondsEnhanced() {
    // Constructor implementation
}

QuantLibBondsEnhanced::~QuantLibBondsEnhanced() {
    // Destructor implementation
}

double QuantLibBondsEnhanced::calculateAmortizingBondPrice(
    const QuantLib::Date& evaluationDate,
    const QuantLib::Date& maturityDate,
    double couponRate,
    QuantLib::Frequency frequency,
    const QuantLib::DayCounter& dayCounter,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    const std::vector<double>& notionals,
    QuantLib::Natural settlementDays) {

    try {
        QuantLib::Settings::instance().evaluationDate() = evaluationDate;

        // Create a schedule
        QuantLib::Calendar calendar = QuantLib::TARGET();
        QuantLib::Date issueDate = evaluationDate - QuantLib::Period(30, QuantLib::Days);

        QuantLib::Schedule schedule(
            issueDate,
            maturityDate,
            QuantLib::Period(frequency),
            calendar,
            QuantLib::Unadjusted,
            QuantLib::Unadjusted,
            QuantLib::DateGeneration::Backward,
            false
        );

        // Set constant coupon rate
        std::vector<double> coupons(notionals.size(), couponRate);

        // Business day convention
        QuantLib::BusinessDayConvention paymentConvention = QuantLib::Following;

        // Create an amortizing bond using the correct constructor signature
        QuantLib::AmortizingFixedRateBond bond(
            settlementDays,
            notionals,
            schedule,
            coupons,
            dayCounter,
            paymentConvention,
            issueDate,               // Pass issueDate as Date, not double
            QuantLib::Period(1, QuantLib::Days),   // Grace period
            calendar,                // Payment calendar
            QuantLib::Following,     // Payment adjustment
            false                    // End of month rule
        );

        // Set pricing engine
        bond.setPricingEngine(
            boost::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::DiscountingBondEngine(discountCurve)
            )
        );

        // Calculate and return the clean price
        return bond.cleanPrice();
    }
    catch (const std::exception& e) {
        std::cerr << "Error in calculateAmortizingBondPrice: " << e.what() << std::endl;
        return 0.0;
    }
}

// Include any other methods from the QuantLibBondsEnhanced class

} // namespace quant