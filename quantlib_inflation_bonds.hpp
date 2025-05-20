// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_inflation_bonds.hpp

#ifndef QUANTLIB_INFLATION_BONDS_HPP
#define QUANTLIB_INFLATION_BONDS_HPP

// Use quotes instead of angle brackets for local includes
#include "ql/quantlib.hpp"
#include <string>
#include <vector>

namespace quant {

class QuantLibInflationBonds {
public:
    // Constructor
    QuantLibInflationBonds(QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed())
        : m_dayCounter(dayCounter) {}

    // Create a zero-coupon inflation indexed bond
    QuantLib::ext::shared_ptr<QuantLib::Bond> createZeroCouponInflationIndexedBond(
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        double baseCPI,
        const QuantLib::Period& observationLag,
        const QuantLib::ext::shared_ptr<QuantLib::ZeroInflationIndex>& inflationIndex,
        QuantLib::Real redemption = 100.0,
        QuantLib::Natural settlementDays = 3);

    // Create a coupon-paying inflation indexed bond
    QuantLib::ext::shared_ptr<QuantLib::Bond> createInflationIndexedBond(
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        double baseCPI,
        const QuantLib::Period& observationLag,
        const QuantLib::ext::shared_ptr<QuantLib::ZeroInflationIndex>& inflationIndex,
        double fixedRate,
        const QuantLib::Period& couponFrequency,
        QuantLib::Real redemption = 100.0,
        QuantLib::Natural settlementDays = 3);

    // Create a zero inflation index from inflation rate curve
    QuantLib::ext::shared_ptr<QuantLib::ZeroInflationIndex> createZeroInflationIndex(
        const std::string& indexName,
        const QuantLib::Handle<QuantLib::ZeroInflationTermStructure>& inflationTS,
        bool interpolated = false);

    // Create a zero inflation curve from market data
    QuantLib::Handle<QuantLib::ZeroInflationTermStructure> createZeroInflationCurve(
        const QuantLib::Date& evaluationDate,
        const std::vector<QuantLib::Date>& dates,
        const std::vector<QuantLib::Rate>& rates,
        double baseRate,
        const QuantLib::Period& observationLag = QuantLib::Period(2, QuantLib::Months));

    // Calculate the inflation-adjusted value at a given date
    double calculateInflationAdjustedValue(
        double notional,
        double baseCPI,
        double currentCPI);

    // Calculate break-even inflation rate between nominal and inflation-indexed bonds
    double calculateBreakEvenInflation(
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& inflationBond,
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& nominalBond,
        double inflationBondPrice,
        double nominalBondPrice);

    // Calculate real yield of an inflation-indexed bond
    double calculateRealYield(
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond,
        double price,
        const QuantLib::Date& settlementDate = QuantLib::Date());

private:
    QuantLib::DayCounter m_dayCounter;
};

} // namespace quant

#endif // QUANTLIB_INFLATION_BONDS_HPP