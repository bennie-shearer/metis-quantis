// File: quantlib_cir_plus_plus_model.hpp

#ifndef QUANTLIB_CIR_PLUS_PLUS_MODEL_H
#define QUANTLIB_CIR_PLUS_PLUS_MODEL_H

#include <ql/quantlib.hpp>

namespace quant {

class QuantLibCIRPlusPlusModel {
public:
    struct CIRPlusPlusParameters {
        double a;       // Mean reversion speed
        double b;       // Long-term mean level
        double sigma;   // Volatility
        double r0;      // Initial short rate
        double h0;      // Initial adjustment term
    };

    // Constructor with default parameters
    QuantLibCIRPlusPlusModel(double a = 0.1,
                            double b = 0.05,
                            double sigma = 0.01,
                            double r0 = 0.03,
                            double h0 = 0.0)
        : m_params({a, b, sigma, r0, h0}) {}

    // Create CIR++ process
    QuantLib::ext::shared_ptr<QuantLib::CoxIngersollRoss> createCIRProcess(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Calibrate model to market data
    void calibrateToMarket(
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
        const std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>>& instruments);

    // Price zero coupon bond using the CIR++ model
    double priceZeroCouponBond(
        const QuantLib::Date& settlementDate,
        const QuantLib::Date& maturityDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Calculate short rate evolution
    std::vector<double> simulateShortRateEvolution(
        const QuantLib::Date& startDate,
        const QuantLib::Date& endDate,
        size_t steps,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Get/set model parameters
    CIRPlusPlusParameters getParameters() const { return m_params; }
    void setParameters(const CIRPlusPlusParameters& params) { m_params = params; }

private:
    CIRPlusPlusParameters m_params;

    // Calculate the deterministic shift function h(t) to match the term structure
    double calculateShift(
        double t,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) const;
};

} // namespace quant

#endif // QUANTLIB_CIR_PLUS_PLUS_MODEL_H