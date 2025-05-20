#ifndef METIS_CALCULATOR_API_HPP
#define METIS_CALCULATOR_API_HPP

#include <string>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/all.hpp>
#include "metis_webserver.hpp"
#include "metis_json.hpp"

// Utility functions - declared as inline to avoid multiple definition errors
inline QuantLib::Date parseDate(const std::string& dateStr) {
    int year, month, day;
    sscanf(dateStr.c_str(), "%d-%d-%d", &year, &month, &day);
    return QuantLib::Date(day, static_cast<QuantLib::Month>(month), year);
}

inline QuantLib::DayCounter parseDayCounter(const std::string& dcStr) {
    if (dcStr == "Actual360")
        return QuantLib::Actual360();
    else if (dcStr == "Actual365Fixed")
        return QuantLib::Actual365Fixed();
    else if (dcStr == "ActualActual")
        return QuantLib::ActualActual(QuantLib::ActualActual::ISDA); // Use ISDA convention
    else
        return QuantLib::Actual360(); // Default
}

inline std::string formatDate(const QuantLib::Date& date) {
    char buffer[11];
    sprintf(buffer, "%04d-%02d-%02d", date.year(), static_cast<int>(date.month()), date.dayOfMonth());
    return std::string(buffer);
}

// API function declarations
void handleApiRequest(const simple_web::Request& request, simple_web::Response& response);
void registerQuantLibCalculatorAPI(simple_web::SimpleWebServer& server);
void registerYieldCurveAPI(simple_web::SimpleWebServer& server);

// Calculator function declarations
simple_json::SimpleJSON calculateBond(const simple_json::SimpleJSON& json);
simple_json::SimpleJSON calculateOption(const simple_json::SimpleJSON& json);
simple_json::SimpleJSON calculateSwap(const simple_json::SimpleJSON& json);
simple_json::SimpleJSON buildYieldCurve(const simple_json::SimpleJSON& json);
simple_json::SimpleJSON calibrateVasicekModel(const simple_json::SimpleJSON& json);
simple_json::SimpleJSON processYieldCurveData(const simple_json::SimpleJSON& json);

#endif // METIS_CALCULATOR_API_HPP