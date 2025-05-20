// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_caps_floors.hpp
#ifndef QUANTLIB_CAPS_FLOORS_H
#define QUANTLIB_CAPS_FLOORS_H

#include <ql/quantlib.hpp>

namespace quant {

    class QuantLibCapsFloors {
    public:
        // Price interest rate cap
        double priceCap(
            const QuantLib::Date& startDate,
            const QuantLib::Date& maturityDate,
            double strike,
            double volatility,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
            QuantLib::Frequency frequency = QuantLib::Quarterly);

        // Price interest rate floor
        double priceFloor(
            const QuantLib::Date& startDate,
            const QuantLib::Date& maturityDate,
            double strike,
            double volatility,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
            QuantLib::Frequency frequency = QuantLib::Quarterly);

        // Price interest rate collar
        double priceCollar(
            const QuantLib::Date& startDate,
            const QuantLib::Date& maturityDate,
            double capStrike,
            double floorStrike,
            double volatility,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
            QuantLib::Frequency frequency = QuantLib::Quarterly);

        // Calculate cap/floor Greeks
        QuantLib::CapFloor::results calculateCapFloorGreeks(
            const QuantLib::Date& startDate,
            const QuantLib::Date& maturityDate,
            double strike,
            double volatility,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
            bool isCap = true,
            QuantLib::Frequency frequency = QuantLib::Quarterly);
    };

} // namespace quant

#endif // QUANTLIB_CAPS_FLOORS_H