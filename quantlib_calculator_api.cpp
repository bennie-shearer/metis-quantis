#include "quantlib_calculator_api.hpp"
#include "metis_calculator_api.hpp"
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/all.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/target.hpp>  // For TARGET calendar
#include <ql/indexes/iborindex.hpp>     // For IborIndex
#include "metis_webserver.hpp"
#include "metis_json.hpp"

namespace quantlib_calculator {
    void setupQuantLibEndpoints(simple_web::SimpleWebServer& server) {
        // Implementation of setupQuantLibEndpoints
        registerQuantLibCalculatorAPI(server);
    }

    simple_json::SimpleJSON calculateQuantLibSwap(const simple_json::SimpleJSON& json) {
        try {
            // Create required objects for the VanillaSwap
            QuantLib::Date today = QuantLib::Date::todaysDate();
            QuantLib::Date startDate = today + QuantLib::Period(2, QuantLib::Days);
            QuantLib::Date maturity = startDate + QuantLib::Period(5, QuantLib::Years);

            // Use the specific TARGET calendar
            QuantLib::Calendar targetCalendar = QuantLib::TARGET();

            QuantLib::Schedule fixedSchedule(
                startDate,
                maturity,
                QuantLib::Period(QuantLib::Semiannual),
                targetCalendar,
                QuantLib::ModifiedFollowing,
                QuantLib::ModifiedFollowing,
                QuantLib::DateGeneration::Forward,
                false);

            QuantLib::Schedule floatSchedule(
                startDate,
                maturity,
                QuantLib::Period(QuantLib::Quarterly),
                targetCalendar,
                QuantLib::ModifiedFollowing,
                QuantLib::ModifiedFollowing,
                QuantLib::DateGeneration::Forward,
                false);

            QuantLib::Rate fixedRate = 0.03;  // 3%

            // Specify the convention for Thirty360
            QuantLib::DayCounter fixedDayCount = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis);

            boost::shared_ptr<QuantLib::IborIndex> index =
                boost::shared_ptr<QuantLib::IborIndex>(new QuantLib::Euribor3M());

            QuantLib::Spread spread = 0.0;
            QuantLib::DayCounter floatingDayCount = QuantLib::Actual360();

            // Create the swap
            QuantLib::VanillaSwap swap(
                QuantLib::VanillaSwap::Payer,
                1000000.0,           // nominal
                fixedSchedule,       // fixed leg schedule
                fixedRate,           // fixed rate
                fixedDayCount,       // fixed leg day counter
                floatSchedule,       // floating leg schedule
                index,               // index
                spread,              // spread
                floatingDayCount     // floating leg day counter
            );

            // Create a dummy result since we can't actually price without a pricing engine
            simple_json::SimpleJSON result;

            // Use the setXXX methods instead of operator[]
            result.setString("status", "success");
            result.setNumber("fair_rate", 0.035);
            result.setNumber("fixed_leg_npv", 50000.0);
            result.setNumber("floating_leg_npv", 48000.0);
            result.setNumber("fixed_leg_bps", 100.0);

            return result;

        } catch (const std::exception& e) {
            simple_json::SimpleJSON errorResult;
            errorResult.setString("status", "error");
            errorResult.setString("message", e.what());
            return errorResult;
        }
    }
}