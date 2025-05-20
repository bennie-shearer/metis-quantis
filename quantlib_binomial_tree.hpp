// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_binomial_tree.hpp

#ifndef QUANT_BOOST_QUANTLIB_BINOMIAL_TREE_HPP
#define QUANT_BOOST_QUANTLIB_BINOMIAL_TREE_HPP

#include <ql/quantlib.hpp>
#include "quantlib_option_calculator.hpp"

namespace quant {

    class QuantLibBinomialTree {
    public:
        QuantLibBinomialTree(int steps);
        ~QuantLibBinomialTree() = default;

        double priceAmericanOption(
            double spotPrice,
            double strikePrice,
            double riskFreeRate,
            double dividendYield,
            double volatility,
            double timeToMaturity,
            QuantLib::Option::Type optionType);

        // Get Greeks for American options
        Greeks getGreeks(
            double spotPrice,
            double strikePrice,
            double riskFreeRate,
            double dividendYield,
            double volatility,
            double timeToMaturity,
            QuantLib::Option::Type optionType,
            QuantLib::Exercise::Type exerciseType);

    private:
        int m_steps;
        QuantLib::DayCounter m_dayCounter;
        QuantLib::Calendar m_calendar;
    };

} // namespace quant

#endif // QUANT_BOOST_QUANTLIB_BINOMIAL_TREE_HPP