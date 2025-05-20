// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_bond_calculator.hpp

#ifndef QUANTLIB_BOND_CALCULATOR_HPP
#define QUANTLIB_BOND_CALCULATOR_HPP

#include "quantlib_includes.hpp"
#include <string>
#include <vector>

namespace quant {

class QuantLibBondCalculator {
public:
    // Constructor
    QuantLibBondCalculator(QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed());

    // Calculate fixed rate bond price and related metrics
    std::pair<double, double> calculateFixedRateBondPrice(
        double faceValue,
        double couponRate,
        double yield,
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Date& settlementDate,
        QuantLib::Frequency frequency = QuantLib::Semiannual,
        const QuantLib::DayCounter& dayCounter = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
        QuantLib::BusinessDayConvention convention = QuantLib::Unadjusted,
        double redemption = 100.0);

    // Calculate zero coupon bond price
    double calculateZeroCouponBondPrice(
        double faceValue,
        double yield,
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Date& settlementDate,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed(),
        QuantLib::BusinessDayConvention convention = QuantLib::Unadjusted,
        double redemption = 100.0);

    // Calculate floating rate bond price
    double calculateFloatingRateBondPrice(
        double faceValue,
        double spread,
        double indexValue,
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Date& settlementDate,
        QuantLib::Frequency frequency = QuantLib::Quarterly,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual360(),
        QuantLib::BusinessDayConvention convention = QuantLib::ModifiedFollowing,
        double redemption = 100.0);

    // Calculate inflation-linked bond price
    double calculateInflationLinkedBondPrice(
        double faceValue,
        double couponRate,
        double indexRatio,
        double realYield,
        const QuantLib::Date& issueDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Date& settlementDate,
        QuantLib::Frequency frequency = QuantLib::Semiannual,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed(),
        QuantLib::BusinessDayConvention convention = QuantLib::Unadjusted,
        double redemption = 100.0);

    // Calculate bond duration and convexity
    std::pair<double, double> calculateBondDurationAndConvexity(
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond,
        double yield,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed());

    // Calculate bond yield to maturity
    double calculateBondYield(
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond,
        double cleanPrice,
        const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed(),
        QuantLib::Compounding compounding = QuantLib::Compounded,
        QuantLib::Frequency frequency = QuantLib::Annual);

    // Generate bond cash flows
    std::vector<std::pair<QuantLib::Date, double>> generateBondCashflows(
        const QuantLib::ext::shared_ptr<QuantLib::Bond>& bond);

private:
    QuantLib::DayCounter m_dayCounter;
    QuantLib::Calendar m_calendar;
};

} // namespace quant

#endif // QUANTLIB_BOND_CALCULATOR_HPP