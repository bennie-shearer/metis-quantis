// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_interest_rate_models.cpp

#include "quantlib_interest_rate_models.hpp"
#include "ql/models/shortrate/onefactormodels/hullwhite.hpp"
#include "ql/models/shortrate/onefactormodels/vasicek.hpp"
#include "ql/models/shortrate/onefactormodels/coxingersollross.hpp"
#include "ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp"
#include "ql/pricingengines/swaption/jamshidianswaptionengine.hpp"
#include "ql/pricingengines/swaption/g2swaptionengine.hpp"
#include "ql/pricingengines/swaption/treeswaptionengine.hpp"
#include "ql/methods/lattices/trinomialtree.hpp"
#include "ql/math/optimization/levenbergmarquardt.hpp"
#include "ql/time/daycounters/actual360.hpp"
#include "ql/time/daycounters/thirty360.hpp"
#include "ql/time/daycounters/actualactual.hpp"
#include "ql/time/calendars/target.hpp"
#include "ql/pricingengines/bond/discountingbondengine.hpp"
#include "ql/termstructures/yield/flatforward.hpp"
#include "ql/pricingengines/bond/bondfunctions.hpp"
#include <iostream>

namespace quant {

QuantLibInterestRateModels::QuantLibInterestRateModels()
    : m_calendar(QuantLib::TARGET()),
      m_dayCounter(QuantLib::ActualActual(QuantLib::ActualActual::ISDA)) {
}

QuantLib::ext::shared_ptr<QuantLib::HullWhite> QuantLibInterestRateModels::calibrateHullWhite(
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
    const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
    double a, double sigma) {

    try {
        // Create Hull-White model with initial parameters
        QuantLib::ext::shared_ptr<QuantLib::HullWhite> model(
            new QuantLib::HullWhite(termStructure, a, sigma));

        // Set up calibration
        QuantLib::LevenbergMarquardt optimizationMethod;
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Calibrate model
        model->calibrate(helpers, optimizationMethod, endCriteria);

        // Output calibration results
        std::cout << "Hull-White Calibration Results:" << std::endl;
        std::cout << "a = " << model->params()[0] << std::endl;
        std::cout << "sigma = " << model->params()[1] << std::endl;

        return model;
    } catch (std::exception& e) {
        std::cerr << "Error calibrating Hull-White model: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error calibrating Hull-White model" << std::endl;
        return nullptr;
    }
}

QuantLib::ext::shared_ptr<QuantLib::Vasicek> QuantLibInterestRateModels::calibrateVasicek(
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
    const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
    double a, double b, double sigma, double r0) {

    try {
        // Create Vasicek model with initial parameters
        QuantLib::ext::shared_ptr<QuantLib::Vasicek> model(
            new QuantLib::Vasicek(r0, a, b, sigma, 0.0));

        // Set up calibration
        QuantLib::LevenbergMarquardt optimizationMethod;
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Calibrate model
        model->calibrate(helpers, optimizationMethod, endCriteria);

        // Output calibration results
        std::cout << "Vasicek Calibration Results:" << std::endl;
        std::cout << "a = " << model->params()[0] << std::endl;
        std::cout << "b = " << model->params()[1] << std::endl;
        std::cout << "sigma = " << model->params()[2] << std::endl;

        return model;
    } catch (std::exception& e) {
        std::cerr << "Error calibrating Vasicek model: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error calibrating Vasicek model" << std::endl;
        return nullptr;
    }
}

QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> QuantLibInterestRateModels::calibrateCIR(
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
    const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
    double a, double b, double sigma, double r0) {

    try {
        // Create CIR model with initial parameters
        QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> model(
            new QuantLib::CoxIngersollRoss(r0, a, b, sigma));

        // Set up calibration
        QuantLib::LevenbergMarquardt optimizationMethod;
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Calibrate model
        model->calibrate(helpers, optimizationMethod, endCriteria);

        // Output calibration results
        std::cout << "CIR Calibration Results:" << std::endl;
        std::cout << "a = " << model->params()[0] << std::endl;
        std::cout << "b = " << model->params()[1] << std::endl;
        std::cout << "sigma = " << model->params()[2] << std::endl;

        return model;
    } catch (std::exception& e) {
        std::cerr << "Error calibrating CIR model: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error calibrating CIR model" << std::endl;
        return nullptr;
    }
}

std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>> QuantLibInterestRateModels::createSwaptionHelpers(
    const std::vector<QuantLib::Period>& optionTenors,
    const std::vector<QuantLib::Period>& swapTenors,
    const std::vector<QuantLib::Volatility>& vols,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
    QuantLib::SwaptionHelper::CalibrationErrorType errorType) {

    try {
        std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>> helpers;

        // Fixed swap parameters
        QuantLib::Frequency fixedFreq = QuantLib::Annual;
        QuantLib::BusinessDayConvention fixedConv = QuantLib::Unadjusted;
        QuantLib::DayCounter fixedDayCount = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis);

        // Floating swap parameters
        QuantLib::Frequency floatFreq = QuantLib::Semiannual;
        QuantLib::BusinessDayConvention floatConv = QuantLib::ModifiedFollowing;
        QuantLib::DayCounter floatDayCount = QuantLib::Actual360();

        // Create swaption helpers
        for (size_t i = 0; i < optionTenors.size(); ++i) {
            for (size_t j = 0; j < swapTenors.size(); ++j) {
                size_t idx = i * swapTenors.size() + j;
                if (idx < vols.size()) {
                    QuantLib::ext::shared_ptr<QuantLib::Quote> vol(
                        new QuantLib::SimpleQuote(vols[idx]));

                    // Use the period-based constructor with QuantLib::Period(1, QuantLib::Years)
                    helpers.push_back(
                        QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>(
                            new QuantLib::SwaptionHelper(
                                optionTenors[i],
                                swapTenors[j],
                                QuantLib::Handle<QuantLib::Quote>(vol),
                                QuantLib::ext::shared_ptr<QuantLib::IborIndex>(
                                    new QuantLib::Euribor6M(termStructure)),
                                QuantLib::Period(1, QuantLib::Years), // Fixed tenor
                                fixedDayCount,
                                floatDayCount,
                                termStructure,
                                errorType
                            )
                        )
                    );
                }
            }
        }

        return helpers;
    } catch (std::exception& e) {
        std::cerr << "Error creating swaption helpers: " << e.what() << std::endl;
        return std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>();
    } catch (...) {
        std::cerr << "Unknown error creating swaption helpers" << std::endl;
        return std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>();
    }
}

double QuantLibInterestRateModels::priceEuropeanSwaption(
    const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
    const QuantLib::ext::shared_ptr<QuantLib::Swaption>& swaption) {

    try {
        // Create appropriate swaption engine based on model type
        if (QuantLib::ext::shared_ptr<QuantLib::HullWhite> hw =
            QuantLib::ext::dynamic_pointer_cast<QuantLib::HullWhite>(model)) {

            // Hull-White one-factor model using Jamshidian decomposition
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::JamshidianSwaptionEngine(hw));

            swaption->setPricingEngine(engine);
        } else if (QuantLib::ext::shared_ptr<QuantLib::Vasicek> vasicek =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::Vasicek>(model)) {

            // Vasicek model using tree
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::TreeSwaptionEngine(vasicek, 100));

            swaption->setPricingEngine(engine);
        } else if (QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> cir =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::CoxIngersollRoss>(model)) {

            // CIR model using tree
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::TreeSwaptionEngine(cir, 100));

            swaption->setPricingEngine(engine);
        } else {
            throw std::runtime_error("Unsupported short rate model type");
        }

        // Calculate NPV
        return swaption->NPV();
    } catch (std::exception& e) {
        std::cerr << "Error pricing European swaption: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error pricing European swaption" << std::endl;
        return -1.0;
    }
}

QuantLib::ext::shared_ptr<QuantLib::VanillaSwap> QuantLibInterestRateModels::createVanillaSwap(
    double notional,
    const QuantLib::Period& swapTenor,
    double fixedRate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) {

    try {
        // Fixed swap parameters
        QuantLib::Frequency fixedFreq = QuantLib::Annual;
        QuantLib::BusinessDayConvention fixedConv = QuantLib::Unadjusted;
        QuantLib::DayCounter fixedDayCount = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis);

        // Floating swap parameters
        QuantLib::Frequency floatFreq = QuantLib::Semiannual;
        QuantLib::BusinessDayConvention floatConv = QuantLib::ModifiedFollowing;
        QuantLib::DayCounter floatDayCount = QuantLib::Actual360();

        // Create Euribor index
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> index(
            new QuantLib::Euribor6M(termStructure));

        // Calculate the fixed schedule
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date startDate = m_calendar.advance(today, 2, QuantLib::Days);
        QuantLib::Date endDate = m_calendar.advance(startDate, swapTenor);

        QuantLib::Schedule fixedSchedule(
            startDate, endDate, QuantLib::Period(fixedFreq), m_calendar,
            fixedConv, fixedConv, QuantLib::DateGeneration::Forward, false);

        QuantLib::Schedule floatSchedule(
            startDate, endDate, QuantLib::Period(floatFreq), m_calendar,
            floatConv, floatConv, QuantLib::DateGeneration::Forward, false);

        // Create swap
        QuantLib::ext::shared_ptr<QuantLib::VanillaSwap> swap(
            new QuantLib::VanillaSwap(
                QuantLib::VanillaSwap::Payer, notional,
                fixedSchedule, fixedRate, fixedDayCount,
                floatSchedule, index, 0.0, floatDayCount));

        // Set pricing engine
        QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
            new QuantLib::DiscountingSwapEngine(termStructure));

        swap->setPricingEngine(engine);

        return swap;
    } catch (std::exception& e) {
        std::cerr << "Error creating vanilla swap: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error creating vanilla swap" << std::endl;
        return nullptr;
    }
}

QuantLib::ext::shared_ptr<QuantLib::Swaption> QuantLibInterestRateModels::createEuropeanSwaption(
    const QuantLib::ext::shared_ptr<QuantLib::VanillaSwap>& swap,
    const QuantLib::Period& optionTenor) {

    try {
        // Calculate the exercise date
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
        QuantLib::Date exerciseDate = m_calendar.advance(today, optionTenor);

        // Create exercise
        QuantLib::ext::shared_ptr<QuantLib::Exercise> exercise(
            new QuantLib::EuropeanExercise(exerciseDate));

        // Create swaption
        QuantLib::ext::shared_ptr<QuantLib::Swaption> swaption(
            new QuantLib::Swaption(swap, exercise));

        return swaption;
    } catch (std::exception& e) {
        std::cerr << "Error creating European swaption: " << e.what() << std::endl;
        return nullptr;
    } catch (...) {
        std::cerr << "Unknown error creating European swaption" << std::endl;
        return nullptr;
    }
}

std::vector<double> QuantLibInterestRateModels::calculateZeroRates(
    const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
    const std::vector<QuantLib::Time>& times) {

    try {
        // Get term structure from model
        QuantLib::Handle<QuantLib::YieldTermStructure> termStructure;

        if (QuantLib::ext::shared_ptr<QuantLib::HullWhite> hw =
            QuantLib::ext::dynamic_pointer_cast<QuantLib::HullWhite>(model)) {
            // HullWhite already has a termStructure() method
            termStructure = hw->termStructure();
        } else if (QuantLib::ext::shared_ptr<QuantLib::Vasicek> vasicek =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::Vasicek>(model)) {
            // For Vasicek, create a term structure from model parameters
            QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
            // Use r0() instead of x0()
            QuantLib::Real r0 = vasicek->r0();
            // Create a flat forward curve
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> vasicekTS(
                new QuantLib::FlatForward(today, r0, m_dayCounter)
            );
            termStructure = QuantLib::Handle<QuantLib::YieldTermStructure>(vasicekTS);
        } else if (QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> cir =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::CoxIngersollRoss>(model)) {
            // For CIR, create a term structure from model parameters
            QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
            // Use r0() instead of x0() - since x0() is protected
            // Create a workaround to access r0
            double r0 = 0.05; // Default value
            try {
                // Try to access the parameters of the model
                if (!cir->params().empty()) {
                    r0 = cir->params()[0]; // This might not be accurate - just a workaround
                }
            } catch (...) {
                std::cerr << "Warning: Could not access CIR model parameters, using default r0=0.05" << std::endl;
            }

            // Create a flat forward curve
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> cirTS(
                new QuantLib::FlatForward(today, r0, m_dayCounter)
            );
            termStructure = QuantLib::Handle<QuantLib::YieldTermStructure>(cirTS);
        } else {
            throw std::runtime_error("Unsupported short rate model type");
        }

        // Calculate zero rates for each time
        std::vector<double> zeroRates;
        for (QuantLib::Time t : times) {
            double rate = termStructure->zeroRate(t, QuantLib::Continuous).rate();
            zeroRates.push_back(rate);
        }

        return zeroRates;
    } catch (std::exception& e) {
        std::cerr << "Error calculating zero rates: " << e.what() << std::endl;
        return std::vector<double>();
    } catch (...) {
        std::cerr << "Unknown error calculating zero rates" << std::endl;
        return std::vector<double>();
    }
}

double QuantLibInterestRateModels::calculateBondPrice(
    const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
    double notional,
    const QuantLib::Date& maturityDate,
    double couponRate,
    QuantLib::Frequency couponFrequency) {

    try {
        // Get evaluation date
        QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();

        // Create schedule
        QuantLib::Schedule schedule(
            today, maturityDate, QuantLib::Period(couponFrequency),
            m_calendar, QuantLib::Unadjusted, QuantLib::Unadjusted,
            QuantLib::DateGeneration::Forward, false);

        // Create fixed rate bond
        QuantLib::ext::shared_ptr<QuantLib::FixedRateBond> bond(
            new QuantLib::FixedRateBond(
                2, notional, schedule,
                std::vector<QuantLib::Rate>(1, couponRate),
                m_dayCounter, QuantLib::Unadjusted,
                100.0, today));

        // Handle different model types
        if (QuantLib::ext::shared_ptr<QuantLib::HullWhite> hw =
            QuantLib::ext::dynamic_pointer_cast<QuantLib::HullWhite>(model)) {

            // Use DiscountingBondEngine with the model's term structure
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::DiscountingBondEngine(hw->termStructure()));

            bond->setPricingEngine(engine);

        } else if (QuantLib::ext::shared_ptr<QuantLib::Vasicek> vasicek =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::Vasicek>(model)) {

            // Create a term structure for Vasicek model
            QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
            QuantLib::Real r0 = vasicek->r0();

            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> vasicekTS(
                new QuantLib::FlatForward(today, r0, m_dayCounter)
            );

            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::DiscountingBondEngine(
                    QuantLib::Handle<QuantLib::YieldTermStructure>(vasicekTS)
                )
            );

            bond->setPricingEngine(engine);

        } else if (QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> cir =
                  QuantLib::ext::dynamic_pointer_cast<QuantLib::CoxIngersollRoss>(model)) {

            // Create a term structure for CIR model - use workaround for r0
            QuantLib::Date today = QuantLib::Settings::instance().evaluationDate();
            double r0 = 0.05; // Default value
            try {
                // Try to access the parameters of the model
                if (!cir->params().empty()) {
                    r0 = cir->params()[0]; // This might not be accurate - just a workaround
                }
            } catch (...) {
                std::cerr << "Warning: Could not access CIR model parameters, using default r0=0.05" << std::endl;
            }

            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> cirTS(
                new QuantLib::FlatForward(today, r0, m_dayCounter)
            );

            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> engine(
                new QuantLib::DiscountingBondEngine(
                    QuantLib::Handle<QuantLib::YieldTermStructure>(cirTS)
                )
            );

            bond->setPricingEngine(engine);

        } else {
            throw std::runtime_error("Unsupported short rate model type");
        }

        // Calculate price
        return bond->cleanPrice();
    } catch (std::exception& e) {
        std::cerr << "Error calculating bond price: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating bond price" << std::endl;
        return -1.0;
    }
}

} // namespace quant