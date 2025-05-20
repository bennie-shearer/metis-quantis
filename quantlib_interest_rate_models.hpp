// Path: C:\Users\bshearer\CLionProjects\quant-boost\quantlib_interest_rate_models.hpp

#ifndef QUANTLIB_INTEREST_RATE_MODELS_HPP
#define QUANTLIB_INTEREST_RATE_MODELS_HPP

// Use quotes instead of angle brackets for local includes
#include "ql/quantlib.hpp"
#include <vector>
#include <string>

namespace quant {

class QuantLibInterestRateModels {
public:
    // Constructor
    QuantLibInterestRateModels();

    // Calibrate Hull-White model
    QuantLib::ext::shared_ptr<QuantLib::HullWhite> calibrateHullWhite(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
        const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
        double a, double sigma);

    // Calibrate Vasicek model
    QuantLib::ext::shared_ptr<QuantLib::Vasicek> calibrateVasicek(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
        const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
        double a, double b, double sigma, double r0);

    // Calibrate CIR model
    QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> calibrateCIR(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
        const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
        double a, double b, double sigma, double r0);

    // Create swaption helpers for calibration
    std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>> createSwaptionHelpers(
        const std::vector<QuantLib::Period>& optionTenors,
        const std::vector<QuantLib::Period>& swapTenors,
        const std::vector<QuantLib::Volatility>& vols,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
        QuantLib::SwaptionHelper::CalibrationErrorType errorType = QuantLib::BlackCalibrationHelper::RelativePriceError);

    // Price a European swaption
    double priceEuropeanSwaption(
        const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
        const QuantLib::ext::shared_ptr<QuantLib::Swaption>& swaption);

    // Create a vanilla swap
    QuantLib::ext::shared_ptr<QuantLib::VanillaSwap> createVanillaSwap(
        double notional,
        const QuantLib::Period& swapTenor,
        double fixedRate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Create a European swaption
    QuantLib::ext::shared_ptr<QuantLib::Swaption> createEuropeanSwaption(
        const QuantLib::ext::shared_ptr<QuantLib::VanillaSwap>& swap,
        const QuantLib::Period& optionTenor);

    // Calculate zero rates from a model
    std::vector<double> calculateZeroRates(
        const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
        const std::vector<QuantLib::Time>& times);

    // Calculate bond price using a model
    double calculateBondPrice(
        const QuantLib::ext::shared_ptr<QuantLib::ShortRateModel>& model,
        double notional,
        const QuantLib::Date& maturityDate,
        double couponRate,
        QuantLib::Frequency couponFrequency);

private:
    QuantLib::Calendar m_calendar;
    QuantLib::DayCounter m_dayCounter;
};

} // namespace quant

#endif // QUANTLIB_INTEREST_RATE_MODELS_HPP