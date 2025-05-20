// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_caps_floors.cpp
#include "quantlib_caps_floors.hpp"
#include <ql/instruments/capfloor.hpp>
#include <ql/pricingengines/capfloor/blackcapfloorengine.hpp>
// Change this include to match your file name
#include <ql/pricingengines/capfloor/analyticcapfloorengine.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/termstructures/volatility/optionlet/constantoptionletvol.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/cashflows/iborcoupon.hpp>  // Add this for IborLeg
#include <iostream>

namespace quant {

double QuantLibCapsFloors::priceCap(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double strike,
    double volatility,
    const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    QuantLib::Frequency frequency) {

    try {
        // Create the schedule
        QuantLib::Schedule schedule(
            startDate,
            maturityDate,
            QuantLib::Period(frequency),
            index->fixingCalendar(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create the cap/floor engine
        QuantLib::Handle<QuantLib::OptionletVolatilityStructure> vol(
            QuantLib::ext::shared_ptr<QuantLib::OptionletVolatilityStructure>(
                new QuantLib::ConstantOptionletVolatility(
                    0,
                    index->fixingCalendar(),
                    QuantLib::ModifiedFollowing,
                    volatility,
                    index->dayCounter())));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::BlackCapFloorEngine(discountCurve, vol));

        // Create a Leg from the schedule
        QuantLib::Leg floatingLeg = QuantLib::IborLeg(schedule, index)
            .withNotionals(1.0)
            .withPaymentDayCounter(index->dayCounter());

        // Create the cap with the modified constructor
        QuantLib::Cap cap(floatingLeg,
                      std::vector<QuantLib::Rate>(1, strike));

        cap.setPricingEngine(engine);

        // Return NPV
        return cap.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error pricing cap: " << e.what() << std::endl;
        return -1.0;
    }
}

double QuantLibCapsFloors::priceFloor(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double strike,
    double volatility,
    const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    QuantLib::Frequency frequency) {

    try {
        // Create the schedule
        QuantLib::Schedule schedule(
            startDate,
            maturityDate,
            QuantLib::Period(frequency),
            index->fixingCalendar(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create the cap/floor engine
        QuantLib::Handle<QuantLib::OptionletVolatilityStructure> vol(
            QuantLib::ext::shared_ptr<QuantLib::OptionletVolatilityStructure>(
                new QuantLib::ConstantOptionletVolatility(
                    0,
                    index->fixingCalendar(),
                    QuantLib::ModifiedFollowing,
                    volatility,
                    index->dayCounter())));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::BlackCapFloorEngine(discountCurve, vol));

        // Create a Leg from the schedule
        QuantLib::Leg floatingLeg = QuantLib::IborLeg(schedule, index)
            .withNotionals(1.0)
            .withPaymentDayCounter(index->dayCounter());

        // Create the floor with the modified constructor
        QuantLib::Floor floor(floatingLeg,
                         std::vector<QuantLib::Rate>(1, strike));

        floor.setPricingEngine(engine);

        // Return NPV
        return floor.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error pricing floor: " << e.what() << std::endl;
        return -1.0;
    }
}

double QuantLibCapsFloors::priceCollar(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double capStrike,
    double floorStrike,
    double volatility,
    const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    QuantLib::Frequency frequency) {

    try {
        // Create the schedule
        QuantLib::Schedule schedule(
            startDate,
            maturityDate,
            QuantLib::Period(frequency),
            index->fixingCalendar(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create the cap/floor engine
        QuantLib::Handle<QuantLib::OptionletVolatilityStructure> vol(
            QuantLib::ext::shared_ptr<QuantLib::OptionletVolatilityStructure>(
                new QuantLib::ConstantOptionletVolatility(
                    0,
                    index->fixingCalendar(),
                    QuantLib::ModifiedFollowing,
                    volatility,
                    index->dayCounter())));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::BlackCapFloorEngine(discountCurve, vol));

        // Create a Leg from the schedule
        QuantLib::Leg floatingLeg = QuantLib::IborLeg(schedule, index)
            .withNotionals(1.0)
            .withPaymentDayCounter(index->dayCounter());

        // Create the collar as a combination of a cap and a floor
        QuantLib::Cap cap(floatingLeg,
                      std::vector<QuantLib::Rate>(1, capStrike));

        QuantLib::Floor floor(floatingLeg,
                         std::vector<QuantLib::Rate>(1, floorStrike));

        cap.setPricingEngine(engine);
        floor.setPricingEngine(engine);

        // Return NPV (long cap, short floor)
        return cap.NPV() - floor.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error pricing collar: " << e.what() << std::endl;
        return -1.0;
    }
}

QuantLib::CapFloor::results QuantLibCapsFloors::calculateCapFloorGreeks(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double strike,
    double volatility,
    const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    bool isCap,
    QuantLib::Frequency frequency) {

    try {
        // Create the schedule
        QuantLib::Schedule schedule(
            startDate,
            maturityDate,
            QuantLib::Period(frequency),
            index->fixingCalendar(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create the cap/floor engine
        QuantLib::Handle<QuantLib::OptionletVolatilityStructure> vol(
            QuantLib::ext::shared_ptr<QuantLib::OptionletVolatilityStructure>(
                new QuantLib::ConstantOptionletVolatility(
                    0,
                    index->fixingCalendar(),
                    QuantLib::ModifiedFollowing,
                    volatility,
                    index->dayCounter())));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::BlackCapFloorEngine(discountCurve, vol));

        // Create a Leg from the schedule
        QuantLib::Leg floatingLeg = QuantLib::IborLeg(schedule, index)
            .withNotionals(1.0)
            .withPaymentDayCounter(index->dayCounter());

        // Create the instrument
        QuantLib::CapFloor::results resultStruct;

        if (isCap) {
            QuantLib::Cap cap(floatingLeg, std::vector<QuantLib::Rate>(1, strike));
            cap.setPricingEngine(engine);

            // Just get the NPV - your QuantLib version doesn't have delta() and vega() methods
            resultStruct.value = cap.NPV();
        } else {
            QuantLib::Floor floor(floatingLeg, std::vector<QuantLib::Rate>(1, strike));
            floor.setPricingEngine(engine);

            resultStruct.value = floor.NPV();
        }

        return resultStruct;
    } catch (std::exception& e) {
        std::cerr << "Error calculating cap/floor greeks: " << e.what() << std::endl;
        return QuantLib::CapFloor::results();
    }
}

} // namespace quant