#ifndef QUANTLIB_SABR_MODEL_H
#define QUANTLIB_SABR_MODEL_H

#include <ql/quantlib.hpp>

namespace quant {

class QuantLibSABRModel {
public:
    // Constructor
    QuantLibSABRModel(double alpha = 0.2,
                     double beta = 0.5,
                     double nu = 0.3,
                     double rho = -0.3)
        : m_alpha(alpha),
          m_beta(beta),
          m_nu(nu),
          m_rho(rho) {}

    // Calculate implied volatility using Hagan's formula
    double impliedVolatility(double forward,
                           double strike,
                           double timeToExpiry) const;

    // Calculate implied volatility at a vector of strikes
    std::vector<double> impliedVolatilityCurve(double forward,
                                            const std::vector<double>& strikes,
                                            double timeToExpiry) const;

    // Calculate option price using SABR model
    double optionPrice(double forward,
                     double strike,
                     double timeToExpiry,
                     double discountFactor,
                     QuantLib::Option::Type optionType = QuantLib::Option::Call) const;

    // Calculate option Greeks using SABR model
    QuantLib::Greeks optionGreeks(double forward,
                               double strike,
                               double timeToExpiry,
                               double discountFactor,
                               QuantLib::Option::Type optionType = QuantLib::Option::Call) const;

    // Calibrate SABR parameters to market volatilities
    void calibrate(double forward,
                 const std::vector<double>& strikes,
                 const std::vector<double>& marketVols,
                 double timeToExpiry);

    // Getters and setters
    double alpha() const { return m_alpha; }
    double beta() const { return m_beta; }
    double nu() const { return m_nu; }
    double rho() const { return m_rho; }

    void setAlpha(double alpha) { m_alpha = alpha; }
    void setBeta(double beta) { m_beta = beta; }
    void setNu(double nu) { m_nu = nu; }
    void setRho(double rho) { m_rho = rho; }

private:
    double m_alpha;  // Initial volatility
    double m_beta;   // Elasticity parameter (0 <= beta <= 1)
    double m_nu;     // Volatility of volatility
    double m_rho;    // Correlation between asset and volatility

    // Helper function for calibration (objective function)
    class SABRError;
};

} // namespace quant

#endif // QUANTLIB_SABR_MODEL_H