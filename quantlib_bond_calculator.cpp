// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_bond_calculator.cpp

#include "quantlib_bond_calculator.hpp"
#include <ql/instruments/bonds/fixedratebond.hpp>
#include <ql/instruments/bonds/zerocouponbond.hpp>
#include <ql/instruments/bonds/floatingratebond.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include <ql/cashflows/couponpricer.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/target.hpp>
#include <iostream>

namespace quant {

QuantLibBondCalculator::QuantLibBondCalculator(QuantLib::DayCounter dayCounter)
    : m_dayCounter(dayCounter),
      m_calendar(QuantLib::TARGET()) {
}

std::pair<double, double> QuantLibBondCalculator::calculateFixedRateBondPrice(
    double faceValue,
    double couponRate,
    double yield,
    const QuantLib::Date& issueDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Date& settlementDate,
    QuantLib::Frequency frequency,
    const QuantLib::DayCounter& dayCounter,
    QuantLib::BusinessDayConvention convention,
    double redemption) {

    try {
        // Set evaluation date to settlement date
        QuantLib::Settings::instance().evaluationDate() = settlementDate;

        // Create schedule
        QuantLib::Schedule schedule(
            issueDate,
            maturityDate,
            QuantLib::Period(frequency),
            m_calendar,
            convention,
            convention,
            QuantLib::DateGeneration::Backward,
            false);

        // Create fixed rate bond
        QuantLib::FixedRateBond bond(
            0,                          // Settlement days
            faceValue,                  // Face value
            schedule,
            std::vector<QuantLib::Rate>(1, couponRate),
            dayCounter,
            convention,
            redemption);

        // Set up pricing engine with flat yield curve
        QuantLib::Handle<QuantLib::YieldTermStructure> discountCurve(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(
                    settlementDate,
                    yield,
                    dayCounter,
                    QuantLib::Compounded,
                    frequency)));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingBondEngine(discountCurve));

        bond.setPricingEngine(engine);

        // Calculate clean and dirty prices
        double cleanPrice = bond.cleanPrice();
        double dirtyPrice = bond.dirtyPrice();

        return {cleanPrice, dirtyPrice};

    } catch (std::exception& e) {
        std::cerr << "Error calculating fixed rate bond price: " << e.what() << std::endl;
        return {-1.0, -1.0};
    } catch (...) {
        std::cerr << "Unknown error calculating fixed rate bond price" << std::endl;
        return {-1.0, -1.0};
    }
}

double QuantLibBondCalculator::calculateZeroCouponBondPrice(
    double faceValue,
    double yield,
    const QuantLib::Date& issueDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Date& settlementDate,
    const QuantLib::DayCounter& dayCounter,
    QuantLib::BusinessDayConvention convention,
    double redemption) {

    try {
        // Set evaluation date to settlement date
        QuantLib::Settings::instance().evaluationDate() = settlementDate;

        // Create zero coupon bond
        QuantLib::ZeroCouponBond bond(
            0,                  // Settlement days
            m_calendar,
            faceValue,
            maturityDate,
            convention,
            redemption,
            issueDate);

        // Set up pricing engine with flat yield curve
        QuantLib::Handle<QuantLib::YieldTermStructure> discountCurve(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(
                    settlementDate,
                    yield,
                    dayCounter,
                    QuantLib::Compounded,
                    QuantLib::Annual)));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingBondEngine(discountCurve));

        bond.setPricingEngine(engine);

        // Return the clean price
        return bond.cleanPrice();

    } catch (std::exception& e) {
        std::cerr << "Error calculating zero coupon bond price: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating zero coupon bond price" << std::endl;
        return -1.0;
    }
}

double QuantLibBondCalculator::calculateFloatingRateBondPrice(
    double faceValue,
    double spread,
    double indexValue,
    const QuantLib::Date& issueDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Date& settlementDate,
    QuantLib::Frequency frequency,
    const QuantLib::DayCounter& dayCounter,
    QuantLib::BusinessDayConvention convention,
    double redemption) {

    try {
        // Set evaluation date to settlement date
        QuantLib::Settings::instance().evaluationDate() = settlementDate;

        // Create schedule
        QuantLib::Schedule schedule(
            issueDate,
            maturityDate,
            QuantLib::Period(frequency),
            m_calendar,
            convention,
            convention,
            QuantLib::DateGeneration::Backward,
            false);

        // Create flat termstructure for the index
        QuantLib::Handle<QuantLib::YieldTermStructure> forecastCurve(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(
                    settlementDate,
                    indexValue,
                    dayCounter,
                    QuantLib::Compounded,
                    frequency)));

        // Create Euribor index
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> index(
            new QuantLib::Euribor3M(forecastCurve));

        // Create floating rate bond
        QuantLib::FloatingRateBond bond(
            0,                      // Settlement days
            faceValue,
            schedule,
            index,
            dayCounter,
            convention,
            2,                      // Fixing days
            std::vector<QuantLib::Real>(1, 1.0),  // Gearings
            std::vector<QuantLib::Spread>(1, spread),  // Spreads
            std::vector<QuantLib::Rate>(),  // Caps
            std::vector<QuantLib::Rate>(),  // Floors
            false,                  // In Arrears
            redemption);

        // Set up pricing engine with flat yield curve (discount curve)
        QuantLib::Handle<QuantLib::YieldTermStructure> discountCurve(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(
                    settlementDate,
                    indexValue + spread,  // Simple assumption for discount rate
                    dayCounter,
                    QuantLib::Compounded,
                    frequency)));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingBondEngine(discountCurve));

        bond.setPricingEngine(engine);

        // Return the clean price
        return bond.cleanPrice();

    } catch (std::exception& e) {
        std::cerr << "Error calculating floating rate bond price: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating floating rate bond price" << std::endl;
        return -1.0;
    }
}

double QuantLibBondCalculator::calculateInflationLinkedBondPrice(
    double faceValue,
    double couponRate,
    double indexRatio,
    double realYield,
    const QuantLib::Date& issueDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Date& settlementDate,
    QuantLib::Frequency frequency,
    const QuantLib::DayCounter& dayCounter,
    QuantLib::BusinessDayConvention convention,
    double redemption) {

    try {
        // This is a simplified implementation since QuantLib's inflation framework is complex
        // We'll simulate an inflation-linked bond using standard cashflows and adjusting for inflation

        // Set evaluation date to settlement date
        QuantLib::Settings::instance().evaluationDate() = settlementDate;

        // Create schedule
        QuantLib::Schedule schedule(
            issueDate,
            maturityDate,
            QuantLib::Period(frequency),
            m_calendar,
            convention,
            convention,
            QuantLib::DateGeneration::Backward,
            false);

        // Create a standard fixed rate bond as a base
        QuantLib::FixedRateBond bond(
            0,                          // Settlement days
            faceValue,                  // Face value
            schedule,
            std::vector<QuantLib::Rate>(1, couponRate),
            dayCounter,
            convention,
            redemption);

        // Set up pricing engine with flat real yield curve
        QuantLib::Handle<QuantLib::YieldTermStructure> discountCurve(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::FlatForward(
                    settlementDate,
                    realYield,
                    dayCounter,
                    QuantLib::Compounded,
                    frequency)));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingBondEngine(discountCurve));

        bond.setPricingEngine(engine);

        // Get the base clean price
        double baseCleanPrice = bond.cleanPrice();

        // Adjust the clean price for inflation
        double inflationAdjustedPrice = baseCleanPrice * indexRatio;

        return inflationAdjustedPrice;

    } catch (std::exception& e) {
        std::cerr << "Error calculating inflation-linked bond price: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating inflation-linked bond price" << std::endl;
        return -1.0;
    }
}

std::pair<double, double> QuantLibBondCalculator::calculateBondDurationAndConvexity(
    const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond,
    double yield,
    const QuantLib::DayCounter& dayCounter) {

    try {
        // Set up yield and evaluation date
        QuantLib::Date settlementDate = QuantLib::Settings::instance().evaluationDate();

        // Calculate Macaulay duration
        double duration = QuantLib::BondFunctions::duration(
            *bond,
            yield,
            dayCounter,
            QuantLib::Compounded,
            QuantLib::Annual,
            QuantLib::Duration::Modified,  // Added Duration type parameter
            settlementDate);

        // Calculate convexity - note convexity() doesn't take a Duration::Type parameter
        double convexity = QuantLib::BondFunctions::convexity(
            *bond,
            yield,
            dayCounter,
            QuantLib::Compounded,
            QuantLib::Annual,
            settlementDate);

        return {duration, convexity};

    } catch (std::exception& e) {
        std::cerr << "Error calculating bond duration and convexity: " << e.what() << std::endl;
        return {-1.0, -1.0};
    } catch (...) {
        std::cerr << "Unknown error calculating bond duration and convexity" << std::endl;
        return {-1.0, -1.0};
    }
}

double QuantLibBondCalculator::calculateBondYield(
    const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond,
    double cleanPrice,
    const QuantLib::DayCounter& dayCounter,
    QuantLib::Compounding compounding,
    QuantLib::Frequency frequency) {

    try {
        // Set up evaluation date
        QuantLib::Date settlementDate = QuantLib::Settings::instance().evaluationDate();

        // Calculate yield
        return QuantLib::BondFunctions::yield(
            *bond,
            cleanPrice,
            dayCounter,
            compounding,
            frequency,
            settlementDate);

    } catch (std::exception& e) {
        std::cerr << "Error calculating bond yield: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating bond yield" << std::endl;
        return -1.0;
    }
}

std::vector<std::pair<QuantLib::Date, double>> QuantLibBondCalculator::generateBondCashflows(
    const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond) {

    try {
        std::vector<std::pair<QuantLib::Date, double>> cashflows;

        // Get the bond's cashflow
        QuantLib::Leg leg = bond->cashflows();

        for (size_t i = 0; i < leg.size(); ++i) {
            QuantLib::Date paymentDate = leg[i]->date();
            double amount = leg[i]->amount();

            if (paymentDate >= QuantLib::Settings::instance().evaluationDate()) {
                cashflows.emplace_back(paymentDate, amount);
            }
        }

        return cashflows;

    } catch (std::exception& e) {
        std::cerr << "Error generating bond cashflows: " << e.what() << std::endl;
        return {};
    } catch (...) {
        std::cerr << "Unknown error generating bond cashflows" << std::endl;
        return {};
    }
}

} // namespace quant