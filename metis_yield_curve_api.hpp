#ifndef METIS_YIELD_CURVE_API_HPP
#define METIS_YIELD_CURVE_API_HPP

#include <string>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/all.hpp>
#include "metis_webserver.hpp"
#include "metis_json.hpp"

namespace yield_curve {
    void setupYieldCurveEndpoints(simple_web::SimpleWebServer& server);
}

#endif // METIS_YIELD_CURVE_API_HPP