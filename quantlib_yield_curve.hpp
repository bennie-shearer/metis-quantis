// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_yield_curve.hpp

#ifndef QUANT_BOOST_QUANTLIB_YIELD_CURVE_HPP
#define QUANT_BOOST_QUANTLIB_YIELD_CURVE_HPP

#include <ql/quantlib.hpp>
#include <string>
#include <vector>

namespace quant {

    class QuantLibYieldCurve {
    public:
        QuantLibYieldCurve();
        ~QuantLibYieldCurve() = default;

        // Build a zero rate curve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> buildZeroCurve(
            const std::vector<QuantLib::Date>& dates,
            const std::vector<QuantLib::Rate>& rates,
            const QuantLib::DayCounter& dayCounter,
            const std::string& interpolationMethod = "Linear");

        // Build a discount curve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> buildDiscountCurve(
            const std::vector<QuantLib::Date>& dates,
            const std::vector<QuantLib::Rate>& discountFactors,
            const QuantLib::DayCounter& dayCounter,
            const std::string& interpolationMethod = "Linear");

        // Build a forward curve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> buildForwardCurve(
            const std::vector<QuantLib::Date>& dates,
            const std::vector<QuantLib::Rate>& forwardRates,
            const QuantLib::DayCounter& dayCounter,
            const std::string& interpolationMethod = "Linear");

    private:
        QuantLib::Calendar m_calendar;
        QuantLib::Natural m_settlementDays;
    };

} // namespace quant

#endif // QUANT_BOOST_QUANTLIB_YIELD_CURVE_HPP