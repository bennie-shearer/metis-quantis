// quantlib_exotic_options.hpp
#pragma once

#include <ql/quantlib.hpp>

namespace quant {

    // Define the enums needed for the option types
    enum class OptionType { Call, Put };
    enum class BarrierType { UpIn, UpOut, DownIn, DownOut };
    enum class AsianAverageType { Arithmetic, Geometric };

    class QuantLibExoticOptions {
    public:
        // Use declarations without default implementation
        QuantLibExoticOptions();
        ~QuantLibExoticOptions();

        // Add declarations for the implemented methods
        double calculateBarrierOptionPrice(
            double spot,
            double strike,
            double barrier,
            double rebate,
            double volatility,
            double riskFreeRate,
            double dividendYield,
            double timeToMaturity,
            OptionType optionType,
            BarrierType barrierType);

        double calculateAsianOptionPrice(
            double spot,
            double strike,
            double volatility,
            double riskFreeRate,
            double dividendYield,
            double timeToMaturity,
            size_t numTimeSteps,
            size_t numPaths,
            OptionType optionType,
            AsianAverageType averageType);
    };

} // namespace quant