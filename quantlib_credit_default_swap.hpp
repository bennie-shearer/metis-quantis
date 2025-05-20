// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_credit_default_swap_enhanced.hpp

#ifndef QUANTLIB_CREDIT_DEFAULT_SWAP_ENHANCED_HPP
#define QUANTLIB_CREDIT_DEFAULT_SWAP_ENHANCED_HPP

#include <ql/quantlib.hpp>
#include <map>
#include <vector>

namespace quant {

class QuantLibCDSEnhanced {
public:
    // Calculate NPV of a Credit Default Swap
    double calculateCdsNPV(
        const QuantLib::Date& valuationDate,
        const QuantLib::Date& maturityDate,
        double spread,                    // in basis points
        double notional,
        double recoveryRate,
        const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve,
        QuantLib::CreditDefaultSwap::PricingModel model = QuantLib::CreditDefaultSwap::Midpoint);

    // Calculate fair spread of a Credit Default Swap
    double calculateCdsFairSpread(
        const QuantLib::Date& valuationDate,
        const QuantLib::Date& maturityDate,
        double notional,
        double recoveryRate,
        const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve,
        QuantLib::CreditDefaultSwap::PricingModel model = QuantLib::CreditDefaultSwap::Midpoint);

    // Bootstrap default probability curve from CDS spreads
    QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure> bootstrapDefaultCurve(
        const QuantLib::Date& valuationDate,
        const std::vector<QuantLib::Period>& tenors,
        const std::vector<double>& spreads,      // in basis points
        double recoveryRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve);

    // Calculate various risk measures for a CDS (NPV, fair spread, expected loss, sensitivities)
    std::map<std::string, double> calculateCdsRiskMeasures(
        const QuantLib::Date& valuationDate,
        const QuantLib::Date& maturityDate,
        double spread,                    // in basis points
        double notional,
        double recoveryRate,
        const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve);
};

} // namespace quant

#endif // QUANTLIB_CREDIT_DEFAULT_SWAP_ENHANCED_HPP