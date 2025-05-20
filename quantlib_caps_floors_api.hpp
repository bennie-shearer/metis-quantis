// File: quantlib_caps_floors_api.hpp
#ifndef QUANTLIB_CAPS_FLOORS_API_HPP
#define QUANTLIB_CAPS_FLOORS_API_HPP

#include "metis_webserver.hpp"
#include "metis_json.hpp"

namespace quant {

    // Cap/Floor API handler
    simple_json::SimpleJSON handleCapFloorCalculation(const std::string& requestBody);

    // Register Cap/Floor API with the web server
    void registerCapFloorAPI(simple_web::SimpleWebServer& server);

} // namespace quant

#endif // QUANTLIB_CAPS_FLOORS_API_HPP