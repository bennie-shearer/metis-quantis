#ifndef QUANTLIB_CALCULATOR_API_HPP
#define QUANTLIB_CALCULATOR_API_HPP

#include <string>
#include <ql/time/date.hpp>
#include "metis_webserver.hpp"
#include "metis_json.hpp"

namespace quantlib_calculator {
    void setupQuantLibEndpoints(simple_web::SimpleWebServer& server);
    simple_json::SimpleJSON calculateQuantLibSwap(const simple_json::SimpleJSON& json);
}

#endif // QUANTLIB_CALCULATOR_API_HPP