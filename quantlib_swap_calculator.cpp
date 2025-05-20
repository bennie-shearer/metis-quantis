// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_swap_calculator.cpp
#include "quantlib_swap_calculator.hpp"
#include <ql/instruments/vanillaswap.hpp>
#include <ql/instruments/makevanillaswap.hpp>
// Comment out the missing header if not needed
// #include <ql/instruments/basisswap.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/indexes/ibor/usdlibor.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/unitedstates.hpp>
#include <iostream>

namespace quant {

double QuantLibSwapCalculator::calculateVanillaSwapNPV(
    double notional,
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double fixedRate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
    QuantLib::Frequency fixedFrequency,
    QuantLib::Frequency floatingFrequency,
    const QuantLib::Period& floatingTenor,
    const QuantLib::DayCounter& fixedDayCounter,
    const QuantLib::DayCounter& floatDayCounter,
    QuantLib::VanillaSwap::Type type) {

    try {
        // Create schedules
        QuantLib::Schedule fixedSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(fixedFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        QuantLib::Schedule floatSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(floatingFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create Euribor index with the appropriate tenor
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> iborIndex(
            new QuantLib::Euribor(floatingTenor, forecastCurve));

        // Create vanilla swap
        QuantLib::VanillaSwap swap(
            type,
            notional,
            fixedSchedule,
            fixedRate,
            fixedDayCounter,
            floatSchedule,
            iborIndex,
            0.0,
            floatDayCounter);

        // Set pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingSwapEngine(discountCurve));

        swap.setPricingEngine(engine);

        // Return NPV
        return swap.NPV();
    } catch (std::exception& e) {
        std::cerr << "Error calculating vanilla swap NPV: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating vanilla swap NPV" << std::endl;
        return -1.0;
    }
}

double QuantLibSwapCalculator::calculateFairFixedRate(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
    QuantLib::Frequency fixedFrequency,
    QuantLib::Frequency floatingFrequency,
    const QuantLib::Period& floatingTenor,
    const QuantLib::DayCounter& fixedDayCounter,
    const QuantLib::DayCounter& floatDayCounter) {

    try {
        // Create Euribor index with the appropriate tenor
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> iborIndex(
            new QuantLib::Euribor(floatingTenor, forecastCurve));

        // Calculate tenor in years for swap
        double years = QuantLib::Actual365Fixed().yearFraction(startDate, maturityDate);
        // Round to nearest integer number of years
        int tenor = static_cast<int>(std::round(years));

        // Create vanilla swap using MakeVanillaSwap
        QuantLib::VanillaSwap swap = QuantLib::MakeVanillaSwap(
            QuantLib::Period(tenor, QuantLib::Years),
            iborIndex,
            0.0)  // Will be replaced with the fair rate
            .withEffectiveDate(startDate)
            .withTerminationDate(maturityDate)
            .withRule(QuantLib::DateGeneration::Forward)
            .withFixedLegTenor(QuantLib::Period(fixedFrequency))
            .withFloatingLegTenor(QuantLib::Period(floatingFrequency))
            .withFixedLegDayCount(fixedDayCounter)
            .withFloatingLegDayCount(floatDayCounter)
            .withDiscountingTermStructure(discountCurve);

        // Return fair rate
        return swap.fairRate();
    } catch (std::exception& e) {
        std::cerr << "Error calculating fair fixed rate: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating fair fixed rate" << std::endl;
        return -1.0;
    }
}

double QuantLibSwapCalculator::calculateFairSpread(
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
    QuantLib::Frequency fixedFrequency,
    QuantLib::Frequency floatingFrequency,
    const QuantLib::Period& floatingTenor,
    const QuantLib::DayCounter& fixedDayCounter,
    const QuantLib::DayCounter& floatDayCounter,
    double fixedRate) {

    try {
        // Create Euribor index with the appropriate tenor
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> iborIndex(
            new QuantLib::Euribor(floatingTenor, forecastCurve));

        // Calculate tenor in years for swap
        double years = QuantLib::Actual365Fixed().yearFraction(startDate, maturityDate);
        // Round to nearest integer number of years
        int tenor = static_cast<int>(std::round(years));

        // Create vanilla swap using MakeVanillaSwap
        QuantLib::VanillaSwap swap = QuantLib::MakeVanillaSwap(
            QuantLib::Period(tenor, QuantLib::Years),
            iborIndex,
            fixedRate)
            .withEffectiveDate(startDate)
            .withTerminationDate(maturityDate)
            .withRule(QuantLib::DateGeneration::Forward)
            .withFixedLegTenor(QuantLib::Period(fixedFrequency))
            .withFloatingLegTenor(QuantLib::Period(floatingFrequency))
            .withFixedLegDayCount(fixedDayCounter)
            .withFloatingLegDayCount(floatDayCounter)
            .withDiscountingTermStructure(discountCurve);

        // Return fair spread
        return swap.fairSpread();
    } catch (std::exception& e) {
        std::cerr << "Error calculating fair spread: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating fair spread" << std::endl;
        return -1.0;
    }
}

SwapResults QuantLibSwapCalculator::calculateSwapGreeks(
    double notional,
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double fixedRate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
    QuantLib::Frequency fixedFrequency,
    QuantLib::Frequency floatingFrequency,
    const QuantLib::Period& floatingTenor,
    const QuantLib::DayCounter& fixedDayCounter,
    const QuantLib::DayCounter& floatDayCounter,
    QuantLib::VanillaSwap::Type type) {

    try {
        // Create schedules
        QuantLib::Schedule fixedSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(fixedFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        QuantLib::Schedule floatSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(floatingFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Create Euribor index with the appropriate tenor
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> iborIndex(
            new QuantLib::Euribor(floatingTenor, forecastCurve));

        // Create vanilla swap
        QuantLib::VanillaSwap swap(
            type,
            notional,
            fixedSchedule,
            fixedRate,
            fixedDayCounter,
            floatSchedule,
            iborIndex,
            0.0,
            floatDayCounter);

        // Set pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingSwapEngine(discountCurve));

        swap.setPricingEngine(engine);

        // Get the swap results
        SwapResults results;
        results.value = swap.NPV();
        results.legBPS[0] = swap.fixedLegBPS();
        results.legBPS[1] = swap.floatingLegBPS();
        results.legNPV[0] = swap.fixedLegNPV();
        results.legNPV[1] = swap.floatingLegNPV();

        return results;
    } catch (std::exception& e) {
        std::cerr << "Error calculating swap Greeks: " << e.what() << std::endl;
        return SwapResults();
    } catch (...) {
        std::cerr << "Unknown error calculating swap Greeks" << std::endl;
        return SwapResults();
    }
}

double QuantLibSwapCalculator::calculateCrossCurrencySwapNPV(
    double domesticNotional,
    double foreignNotional,
    const QuantLib::Date& startDate,
    const QuantLib::Date& maturityDate,
    double domesticRate,
    double foreignRate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& domesticCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& foreignCurve,
    const QuantLib::Handle<QuantLib::Quote>& fxRate,
    QuantLib::Frequency domesticFrequency,
    QuantLib::Frequency foreignFrequency,
    const QuantLib::DayCounter& domesticDayCounter,
    const QuantLib::DayCounter& foreignDayCounter) {

    try {
        // Create schedules
        QuantLib::Schedule domesticSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(domesticFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        QuantLib::Schedule foreignSchedule(
            startDate,
            maturityDate,
            QuantLib::Period(foreignFrequency),
            QuantLib::TARGET(),
            QuantLib::ModifiedFollowing,
            QuantLib::ModifiedFollowing,
            QuantLib::DateGeneration::Forward,
            false);

        // Calculate domestic and foreign leg NPVs separately

        // Domestic leg
        QuantLib::Leg domesticLeg = QuantLib::FixedRateLeg(domesticSchedule)
            .withNotionals(domesticNotional)
            .withCouponRates(domesticRate, domesticDayCounter)
            .withPaymentAdjustment(QuantLib::ModifiedFollowing);

        // Get the underlying YieldTermStructure reference using currentLink()
        const QuantLib::YieldTermStructure& domesticTS = *(domesticCurve.currentLink());

        // Calculate NPV of the domestic leg
        double domesticLegNPV = QuantLib::CashFlows::npv(
            domesticLeg,
            domesticTS,
            false,             // includeSettlementDateFlows
            domesticCurve->referenceDate(),
            domesticCurve->referenceDate());

        // Foreign leg
        QuantLib::Leg foreignLeg = QuantLib::FixedRateLeg(foreignSchedule)
            .withNotionals(foreignNotional)
            .withCouponRates(foreignRate, foreignDayCounter)
            .withPaymentAdjustment(QuantLib::ModifiedFollowing);

        // Get the underlying YieldTermStructure reference using currentLink()
        const QuantLib::YieldTermStructure& foreignTS = *(foreignCurve.currentLink());

        // Calculate NPV of the foreign leg
        double foreignLegNPV = QuantLib::CashFlows::npv(
            foreignLeg,
            foreignTS,
            false,             // includeSettlementDateFlows
            foreignCurve->referenceDate(),
            foreignCurve->referenceDate());

        // Convert foreign NPV to domestic currency and calculate swap NPV
        double swapNPV = domesticLegNPV - foreignLegNPV * fxRate->value();

        return swapNPV;
    } catch (std::exception& e) {
        std::cerr << "Error calculating cross currency swap NPV: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating cross currency swap NPV" << std::endl;
        return -1.0;
    }
}

} // namespace quant