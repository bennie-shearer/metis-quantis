// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_yield_curve.cpp

#include "quantlib_yield_curve.hpp"
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/termstructures/yield/forwardcurve.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <iostream>

namespace quant {

QuantLibYieldCurve::QuantLibYieldCurve()
    : m_calendar(QuantLib::TARGET()),
      m_settlementDays(2) {
}

QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> QuantLibYieldCurve::buildZeroCurve(
    const std::vector<QuantLib::Date>& dates,
    const std::vector<QuantLib::Rate>& rates,
    const QuantLib::DayCounter& dayCounter,
    const std::string& interpolationMethod) {

    try {
        // Ensure dates and rates have same size
        if (dates.size() != rates.size()) {
            throw std::runtime_error("Dates and rates must have the same size");
        }

        if (dates.empty()) {
            throw std::runtime_error("Dates vector cannot be empty");
        }

        // Create a simple ZeroCurve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> curve(
            new QuantLib::ZeroCurve(dates, rates, dayCounter, m_calendar)
        );

        return curve;
    }
    catch (std::exception& e) {
        std::cerr << "Error building zero curve: " << e.what() << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
    catch (...) {
        std::cerr << "Unknown error building zero curve" << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
}

QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> QuantLibYieldCurve::buildDiscountCurve(
    const std::vector<QuantLib::Date>& dates,
    const std::vector<QuantLib::Rate>& discountFactors,
    const QuantLib::DayCounter& dayCounter,
    const std::string& interpolationMethod) {

    try {
        // Ensure dates and discount factors have same size
        if (dates.size() != discountFactors.size()) {
            throw std::runtime_error("Dates and discount factors must have the same size");
        }

        if (dates.empty()) {
            throw std::runtime_error("Dates vector cannot be empty");
        }

        // Create a simple DiscountCurve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> curve(
            new QuantLib::DiscountCurve(dates, discountFactors, dayCounter, m_calendar)
        );

        return curve;
    }
    catch (std::exception& e) {
        std::cerr << "Error building discount curve: " << e.what() << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
    catch (...) {
        std::cerr << "Unknown error building discount curve" << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
}

QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> QuantLibYieldCurve::buildForwardCurve(
    const std::vector<QuantLib::Date>& dates,
    const std::vector<QuantLib::Rate>& forwardRates,
    const QuantLib::DayCounter& dayCounter,
    const std::string& interpolationMethod) {

    try {
        // Ensure dates and forward rates have same size
        if (dates.size() != forwardRates.size()) {
            throw std::runtime_error("Dates and forward rates must have the same size");
        }

        if (dates.empty()) {
            throw std::runtime_error("Dates vector cannot be empty");
        }

        // Create a simple ForwardCurve
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> curve(
            new QuantLib::ForwardCurve(dates, forwardRates, dayCounter, m_calendar)
        );

        return curve;
    }
    catch (std::exception& e) {
        std::cerr << "Error building forward curve: " << e.what() << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
    catch (...) {
        std::cerr << "Unknown error building forward curve" << std::endl;
        return QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>();
    }
}

} // namespace quant