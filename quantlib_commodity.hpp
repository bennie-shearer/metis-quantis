// C:/Users/bshearer/CLionProjects/quant-boost/quantlib_commodity.hpp

#ifndef QUANT_BOOST_QUANTLIB_COMMODITY_HPP
#define QUANT_BOOST_QUANTLIB_COMMODITY_HPP

#include <ql/instrument.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/instruments/payoffs.hpp>  // Correct path for payoffs in QuantLib 1.30
#include <ql/exercise.hpp>

namespace quant {

    // This class will wrap the QuantLib functionality for commodities
    class QuantLibCommodity {
    public:
        QuantLibCommodity();
        ~QuantLibCommodity();

        // Forward pricing
        double calculateCommodityForwardPrice(
            double spotPrice,
            double contractValue,
            const QuantLib::Date& evaluationDate,
            const QuantLib::Date& maturityDate,
            double riskFreeRate,
            double storageRate,
            double convenienceYield,
            const QuantLib::DayCounter& dayCounter
        );

        // Option pricing
        double calculateCommodityOptionPrice(
            double spotPrice,
            double strikePrice,
            const QuantLib::Date& evaluationDate,
            const QuantLib::Date& maturityDate,
            double riskFreeRate,
            double volatility,
            const QuantLib::DayCounter& dayCounter,
            bool isCall
        );

    private:
        // Internal implementation details
    };

} // namespace quant

#endif // QUANT_BOOST_QUANTLIB_COMMODITY_HPP