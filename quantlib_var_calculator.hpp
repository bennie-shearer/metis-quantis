// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_var_calculator.hpp
#ifndef QUANTLIB_VAR_CALCULATOR_H
#define QUANTLIB_VAR_CALCULATOR_H

#include <ql/quantlib.hpp>
#include "metis_json.hpp"
#include <vector>
#include <string>

namespace quant {

    // Main VaR calculation function
    simple_json::SimpleJSON calculateVaR(const simple_json::SimpleJSON& params);

} // namespace quant

#endif // QUANTLIB_VAR_CALCULATOR_H