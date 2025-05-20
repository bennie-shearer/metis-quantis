// quantlib_includes.hpp
// A central header for QuantLib includes to replace ql/quantlib.hpp

#ifndef QUANTLIB_INCLUDES_HPP
#define QUANTLIB_INCLUDES_HPP

// Core QuantLib headers
#include <ql/qldefines.hpp>
#include <ql/settings.hpp>
#include <ql/utilities/null.hpp>
#include <ql/shared_ptr.hpp>
#include <ql/handle.hpp>
#include <ql/errors.hpp>

// Time-related headers
#include <ql/time/date.hpp>
#include <ql/time/period.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/unitedstates.hpp>
#include <ql/time/calendars/unitedkingdom.hpp>
#include <ql/time/calendars/japan.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/frequency.hpp>
#include <ql/time/businessdayconvention.hpp>

// Interest rate related headers
#include <ql/interestrate.hpp>
#include <ql/compounding.hpp>
#include <ql/cashflow.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/coupon.hpp>

// Instrument and pricing headers
#include <ql/instrument.hpp>
#include <ql/pricingengine.hpp>
#include <ql/termstructure.hpp>
#include <ql/instruments/bond.hpp>
#include <ql/instruments/bonds/fixedratebond.hpp>
#include <ql/instruments/bonds/floatingratebond.hpp>
#include <ql/instruments/bonds/zerocouponbond.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include <ql/pricingengines/bond/bondfunctions.hpp>

// Term structures
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

// Math and utilities
#include <ql/math/matrix.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/indexes/ibor/usdlibor.hpp>

// Models and methods
#include <ql/models/shortrate/onefactormodel.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/methods/montecarlo/multipathgenerator.hpp>

// Add any other headers your specific project needs here

#endif // QUANTLIB_INCLUDES_HPP