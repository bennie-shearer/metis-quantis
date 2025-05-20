// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_credit_default_swap_enhanced.cpp

#include "quantlib_credit_default_swap.hpp"
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/pricingengines/credit/integralcdsengine.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/calendars/weekendsonly.hpp>
#include <ql/quotes/simplequote.hpp>
#include <iostream>

namespace quant {

double QuantLibCDSEnhanced::calculateCdsNPV(
    const QuantLib::Date& valuationDate,
    const QuantLib::Date& maturityDate,
    double spread,
    double notional,
    double recoveryRate,
    const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve,
    QuantLib::CreditDefaultSwap::PricingModel model) {

    try {
        // Set evaluation date
        QuantLib::Settings::instance().evaluationDate() = valuationDate;

        // Create CDS schedule
        QuantLib::Schedule schedule(
            valuationDate,
            maturityDate,
            QuantLib::Period(QuantLib::Quarterly),
            QuantLib::WeekendsOnly(),
            QuantLib::Following,
            QuantLib::Unadjusted,
            QuantLib::DateGeneration::CDS,
            false);

        // Create CDS
        QuantLib::CreditDefaultSwap cds(
            QuantLib::Protection::Buyer,  // Protection buyer
            notional,
            spread / 10000.0,             // Convert from bps to decimal
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,                        // Settles accrual
            true,                        // Pays at default time
            valuationDate);              // Protection start date

        // Create pricing engine based on model choice
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine;

        if (model == QuantLib::CreditDefaultSwap::Midpoint) {
            engine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::MidPointCdsEngine(defaultCurve, recoveryRate, yieldCurve));
        } else {
            engine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::IntegralCdsEngine(
                    QuantLib::Period(1, QuantLib::Days),  // Integration step
                    defaultCurve,
                    recoveryRate,
                    yieldCurve));
        }

        cds.setPricingEngine(engine);

        // Calculate NPV
        return cds.NPV();
    }
    catch (std::exception& e) {
        std::cerr << "Error calculating CDS NPV: " << e.what() << std::endl;
        return -1.0;
    }
}

double QuantLibCDSEnhanced::calculateCdsFairSpread(
    const QuantLib::Date& valuationDate,
    const QuantLib::Date& maturityDate,
    double notional,
    double recoveryRate,
    const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve,
    QuantLib::CreditDefaultSwap::PricingModel model) {

    try {
        // Set evaluation date
        QuantLib::Settings::instance().evaluationDate() = valuationDate;

        // Create CDS schedule
        QuantLib::Schedule schedule(
            valuationDate,
            maturityDate,
            QuantLib::Period(QuantLib::Quarterly),
            QuantLib::WeekendsOnly(),
            QuantLib::Following,
            QuantLib::Unadjusted,
            QuantLib::DateGeneration::CDS,
            false);

        // Create CDS with dummy spread (will be replaced with fair spread)
        QuantLib::CreditDefaultSwap cds(
            QuantLib::Protection::Buyer,  // Protection buyer
            notional,
            0.01,                        // Dummy spread
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,                        // Settles accrual
            true,                        // Pays at default time
            valuationDate);              // Protection start date

        // Create pricing engine based on model choice
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine;

        if (model == QuantLib::CreditDefaultSwap::Midpoint) {
            engine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::MidPointCdsEngine(defaultCurve, recoveryRate, yieldCurve));
        } else {
            engine = QuantLib::ext::shared_ptr<QuantLib::PricingEngine>(
                new QuantLib::IntegralCdsEngine(
                    QuantLib::Period(1, QuantLib::Days),  // Integration step
                    defaultCurve,
                    recoveryRate,
                    yieldCurve));
        }

        cds.setPricingEngine(engine);

        // Calculate fair spread (converted to basis points)
        return cds.fairSpread() * 10000.0;
    }
    catch (std::exception& e) {
        std::cerr << "Error calculating CDS fair spread: " << e.what() << std::endl;
        return -1.0;
    }
}

QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure> QuantLibCDSEnhanced::bootstrapDefaultCurve(
    const QuantLib::Date& valuationDate,
    const std::vector<QuantLib::Period>& tenors,
    const std::vector<double>& spreads,
    double recoveryRate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve) {

    try {
        QuantLib::Settings::instance().evaluationDate() = valuationDate;

        // Create CDS helpers for bootstrapping
        std::vector<QuantLib::ext::shared_ptr<QuantLib::DefaultProbabilityHelper>> instruments;

        for (size_t i = 0; i < tenors.size(); ++i) {
            QuantLib::ext::shared_ptr<QuantLib::Quote> spread(
                new QuantLib::SimpleQuote(spreads[i] / 10000.0));  // Convert from bps to decimal

            instruments.push_back(
                QuantLib::ext::shared_ptr<QuantLib::SpreadCdsHelper>(
                    new QuantLib::SpreadCdsHelper(
                        QuantLib::Handle<QuantLib::Quote>(spread),
                        tenors[i],
                        1,                        // Settlement days
                        QuantLib::WeekendsOnly(),
                        QuantLib::Quarterly,
                        QuantLib::Following,
                        QuantLib::DateGeneration::CDS,
                        QuantLib::Actual360(),
                        recoveryRate,
                        yieldCurve)));
        }

        // Bootstrap default curve
        QuantLib::ext::shared_ptr<QuantLib::DefaultProbabilityTermStructure> defaultCurve(
            new QuantLib::PiecewiseDefaultCurve<QuantLib::SurvivalProbability, QuantLib::Linear>(
                valuationDate, instruments, QuantLib::Actual365Fixed()));

        defaultCurve->enableExtrapolation();

        return QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>(defaultCurve);
    }
    catch (std::exception& e) {
        std::cerr << "Error bootstrapping default curve: " << e.what() << std::endl;

        // Return a flat hazard rate as fallback
        QuantLib::ext::shared_ptr<QuantLib::DefaultProbabilityTermStructure> flatCurve(
            new QuantLib::FlatHazardRate(
                valuationDate,
                0.01,             // Hazard rate of 1%
                QuantLib::Actual365Fixed()));

        return QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>(flatCurve);
    }
}

std::map<std::string, double> QuantLibCDSEnhanced::calculateCdsRiskMeasures(
    const QuantLib::Date& valuationDate,
    const QuantLib::Date& maturityDate,
    double spread,
    double notional,
    double recoveryRate,
    const QuantLib::Handle<QuantLib::DefaultProbabilityTermStructure>& defaultCurve,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldCurve) {

    try {
        // Set evaluation date
        QuantLib::Settings::instance().evaluationDate() = valuationDate;

        // Create CDS schedule
        QuantLib::Schedule schedule(
            valuationDate,
            maturityDate,
            QuantLib::Period(QuantLib::Quarterly),
            QuantLib::WeekendsOnly(),
            QuantLib::Following,
            QuantLib::Unadjusted,
            QuantLib::DateGeneration::CDS,
            false);

        // Create CDS
        QuantLib::CreditDefaultSwap cds(
            QuantLib::Protection::Buyer,  // Protection buyer
            notional,
            spread / 10000.0,             // Convert from bps to decimal
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,                        // Settles accrual
            true,                        // Pays at default time
            valuationDate);              // Protection start date

        // Create pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::MidPointCdsEngine(defaultCurve, recoveryRate, yieldCurve));

        cds.setPricingEngine(engine);

        // Calculate risk measures
        std::map<std::string, double> results;

        // NPV
        results["NPV"] = cds.NPV();

        // Fair spread
        results["fairSpread"] = cds.fairSpread() * 10000.0;  // Convert to bps

        // Expected loss
        double expectedLoss = 0.0;
        for (QuantLib::Date date = valuationDate; date <= maturityDate; date = date + 1) {
            double dt = defaultCurve->defaultProbability(date, date + 1);
            double df = yieldCurve->discount(date);
            expectedLoss += dt * (1.0 - recoveryRate) * notional * df;
        }
        results["expectedLoss"] = expectedLoss;

        // Credit DV01 (spread sensitivity)
        double originalNPV = results["NPV"];

        // Recalculate with slightly higher spread
        QuantLib::CreditDefaultSwap cdsUp(
            QuantLib::Protection::Buyer,  // Protection buyer
            notional,
            (spread + 1.0) / 10000.0,     // 1 bp higher
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,                        // Settles accrual
            true,                        // Pays at default time
            valuationDate);              // Protection start date

        cdsUp.setPricingEngine(engine);
        double upNPV = cdsUp.NPV();

        results["creditDV01"] = upNPV - originalNPV;  // Change in NPV for 1 bp increase

        // Recovery01 (recovery sensitivity)
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engineUp(
            new QuantLib::MidPointCdsEngine(
                defaultCurve, recoveryRate + 0.01, yieldCurve));  // 1% higher recovery

        QuantLib::CreditDefaultSwap cdsRecUp(
            QuantLib::Protection::Buyer,
            notional,
            spread / 10000.0,
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,
            true,
            valuationDate);

        cdsRecUp.setPricingEngine(engineUp);
        results["recovery01"] = cdsRecUp.NPV() - originalNPV;  // Change in NPV for 1% increase in recovery

        // IR01 (interest rate sensitivity)
        QuantLib::ext::shared_ptr<QuantLib::SimpleQuote> rateQuote(
            new QuantLib::SimpleQuote(0.0));

        QuantLib::Handle<QuantLib::Quote> rateHandle(rateQuote);

        QuantLib::Handle<QuantLib::YieldTermStructure> yieldCurveUp(
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>(
                new QuantLib::ZeroSpreadedTermStructure(
                    yieldCurve, rateHandle)));

        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engineIR(
            new QuantLib::MidPointCdsEngine(
                defaultCurve, recoveryRate, yieldCurveUp));

        QuantLib::CreditDefaultSwap cdsIR(
            QuantLib::Protection::Buyer,
            notional,
            spread / 10000.0,
            schedule,
            QuantLib::Following,
            QuantLib::Actual360(),
            true,
            true,
            valuationDate);

        cdsIR.setPricingEngine(engineIR);
        double baseNPV = cdsIR.NPV();

        rateQuote->setValue(0.0001);  // 1 bp increase

        results["IR01"] = cdsIR.NPV() - baseNPV;  // Change in NPV for 1 bp increase in rates

        return results;
    }
    catch (std::exception& e) {
        std::cerr << "Error calculating CDS risk measures: " << e.what() << std::endl;
        return std::map<std::string, double>();
    }
}

} // namespace quant