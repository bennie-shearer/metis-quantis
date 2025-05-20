// C:/Users/bshearer/CLionProjects/quant-boost/quantlib_bonds.hpp

#ifndef QUANT_BOOST_QUANTLIB_BONDS_HPP
#define QUANT_BOOST_QUANTLIB_BONDS_HPP

#include <ql/quantlib.hpp>
#include <vector>

namespace quant {

    class QuantLibBondsEnhanced {
    public:
        QuantLibBondsEnhanced();
        ~QuantLibBondsEnhanced();

        double calculateAmortizingBondPrice(
            const QuantLib::Date& evaluationDate,
            const QuantLib::Date& maturityDate,
            double couponRate,
            QuantLib::Frequency frequency,
            const QuantLib::DayCounter& dayCounter,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve,
            const std::vector<double>& notionals,
            QuantLib::Natural settlementDays = 3
        );

        // Other methods...
    };

} // namespace quant

#endif // QUANT_BOOST_QUANTLIB_BONDS_HPP