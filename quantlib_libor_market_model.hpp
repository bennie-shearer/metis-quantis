// File: quantlib_libor_market_model.hpp

#ifndef QUANTLIB_LIBOR_MARKET_MODEL_H
#define QUANTLIB_LIBOR_MARKET_MODEL_H

#include <ql/quantlib.hpp>
#include <vector>

namespace quant {

class QuantLibLiborMarketModel {
public:
    // Constructor
    QuantLibLiborMarketModel(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure = QuantLib::Handle<QuantLib::YieldTermStructure>());

    // Set up the LIBOR Market Model
    void setupLMM(
        const std::vector<QuantLib::Time>& rateTimes,
        const std::vector<QuantLib::Real>& initialRates,
        const std::vector<QuantLib::Real>& volatilities,
        const std::vector<QuantLib::Real>& correlations);

    // Calibrate the model to caplet or swaption volatilities
    void calibrateToCaplets(
        const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& helpers,
        const QuantLib::Real lambda = 0.0,
        const QuantLib::Size maxIterations = 1000);

    // Price a caplet using the LMM
    QuantLib::Real priceCaplet(
        QuantLib::Time start,
        QuantLib::Time end,
        QuantLib::Rate strike,
        const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& index) const;

    // Price a swaption using the LMM
    QuantLib::Real priceSwaption(
        const QuantLib::ext::shared_ptr<QuantLib::Swaption>& swaption) const;

    // Simulate forward LIBOR rates
    std::vector<std::vector<QuantLib::Real>> simulateLiborRates(
        QuantLib::Size numberOfPaths,
        QuantLib::Size seed = 42) const;

    // Get the LMM model instance
    QuantLib::ext::shared_ptr<QuantLib::LmLinearExponentialVolModel> getModel() const { return m_model; }

    // Get the calibrated volatility parameters
    std::vector<QuantLib::Real> getVolatilities() const;

    // Get the calibrated correlation parameters
    std::vector<QuantLib::Real> getCorrelations() const;

private:
    QuantLib::Handle<QuantLib::YieldTermStructure> m_termStructure;
    QuantLib::ext::shared_ptr<QuantLib::LmLinearExponentialVolModel> m_model;
    QuantLib::ext::shared_ptr<QuantLib::LiborForwardModelProcess> m_process;
    std::vector<QuantLib::Time> m_rateTimes;
};

} // namespace quant

#endif // QUANTLIB_LIBOR_MARKET_MODEL_H