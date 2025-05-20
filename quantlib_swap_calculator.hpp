// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_swap_calculator.hpp
#ifndef QUANTLIB_SWAP_CALCULATOR_H
#define QUANTLIB_SWAP_CALCULATOR_H

#include <ql/quantlib.hpp>

namespace quant {

// Define a custom swap results structure since QuantLib doesn't provide one
struct SwapResults {
    double value = 0.0;
    std::array<double, 2> legBPS = {0.0, 0.0};
    std::array<double, 2> legNPV = {0.0, 0.0};
};

class QuantLibSwapCalculator {
public:
    // Calculate vanilla interest rate swap NPV
    double calculateVanillaSwapNPV(
        double notional,
        const QuantLib::Date& startDate,
        const QuantLib::Date& maturityDate,
        double fixedRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
        QuantLib::Frequency fixedFrequency = QuantLib::Annual,
        QuantLib::Frequency floatingFrequency = QuantLib::Quarterly,
        const QuantLib::Period& floatingTenor = QuantLib::Period(3, QuantLib::Months),
        const QuantLib::DayCounter& fixedDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        const QuantLib::DayCounter& floatDayCounter = QuantLib::Actual360(),
        QuantLib::VanillaSwap::Type type = QuantLib::VanillaSwap::Payer);

    // Calculate fair fixed rate for an interest rate swap
    double calculateFairFixedRate(
        const QuantLib::Date& startDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
        QuantLib::Frequency fixedFrequency = QuantLib::Annual,
        QuantLib::Frequency floatingFrequency = QuantLib::Quarterly,
        const QuantLib::Period& floatingTenor = QuantLib::Period(3, QuantLib::Months),
        const QuantLib::DayCounter& fixedDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        const QuantLib::DayCounter& floatDayCounter = QuantLib::Actual360());

    // Calculate fair spread for an interest rate swap
    double calculateFairSpread(
        const QuantLib::Date& startDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
        QuantLib::Frequency fixedFrequency = QuantLib::Annual,
        QuantLib::Frequency floatingFrequency = QuantLib::Quarterly,
        const QuantLib::Period& floatingTenor = QuantLib::Period(3, QuantLib::Months),
        const QuantLib::DayCounter& fixedDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        const QuantLib::DayCounter& floatDayCounter = QuantLib::Actual360(),
        double fixedRate = 0.05);

    // Calculate swap Greeks (Value, BPS, NPV for each leg)
    SwapResults calculateSwapGreeks(
        double notional,
        const QuantLib::Date& startDate,
        const QuantLib::Date& maturityDate,
        double fixedRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& forecastCurve,
        QuantLib::Frequency fixedFrequency = QuantLib::Annual,
        QuantLib::Frequency floatingFrequency = QuantLib::Quarterly,
        const QuantLib::Period& floatingTenor = QuantLib::Period(3, QuantLib::Months),
        const QuantLib::DayCounter& fixedDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        const QuantLib::DayCounter& floatDayCounter = QuantLib::Actual360(),
        QuantLib::VanillaSwap::Type type = QuantLib::VanillaSwap::Payer);

    // Calculate cross-currency swap NPV
    double calculateCrossCurrencySwapNPV(
        double domesticNotional,
        double foreignNotional,
        const QuantLib::Date& startDate,
        const QuantLib::Date& maturityDate,
        double domesticRate,
        double foreignRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& domesticCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& foreignCurve,
        const QuantLib::Handle<QuantLib::Quote>& fxRate,
        QuantLib::Frequency domesticFrequency = QuantLib::Annual,
        QuantLib::Frequency foreignFrequency = QuantLib::Annual,
        const QuantLib::DayCounter& domesticDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        const QuantLib::DayCounter& foreignDayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis));
};

} // namespace quant

#endif // QUANTLIB_SWAP_CALCULATOR_H