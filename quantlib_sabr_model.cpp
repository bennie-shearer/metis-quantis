// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_sabr_model.cpp
#include "quantlib_sabr_model.hpp"
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

namespace quant {

// Helper class for calibration (objective function)
class QuantLibSABRModel::SABRError : public QuantLib::CostFunction {
public:
    SABRError(double forward,
             const std::vector<double>& strikes,
             const std::vector<double>& marketVols,
             double timeToExpiry,
             bool fixBeta)
        : m_forward(forward),
          m_strikes(strikes),
          m_marketVols(marketVols),
          m_timeToExpiry(timeToExpiry),
          m_fixBeta(fixBeta) {}

    QuantLib::Real value(const QuantLib::Array& params) const override {
        double alpha = params[0];
        double beta = m_fixBeta ? 0.5 : params[1]; // Can fix beta to a common value
        double nu = m_fixBeta ? params[1] : params[2];
        double rho = m_fixBeta ? params[2] : params[3];

        // Ensure parameters are in valid range
        if (alpha <= 0.0 || nu <= 0.0 || rho <= -1.0 || rho >= 1.0 || beta < 0.0 || beta > 1.0) {
            return 1.0e10; // Return a large number instead of QL_MAX_REAL
        }

        QuantLibSABRModel sabr(alpha, beta, nu, rho);

        double error = 0.0;
        for (size_t i = 0; i < m_strikes.size(); ++i) {
            double modelVol = sabr.impliedVolatility(m_forward, m_strikes[i], m_timeToExpiry);
            double marketVol = m_marketVols[i];

            // Use relative error for better fit
            double diff = (modelVol - marketVol) / marketVol;
            error += diff * diff;
        }

        return error;
    }

    // Implement the pure virtual values() method required by CostFunction
    QuantLib::Array values(const QuantLib::Array& params) const override {
        QuantLib::Array result(1);
        result[0] = value(params);
        return result;
    }

private:
    double m_forward;
    std::vector<double> m_strikes;
    std::vector<double> m_marketVols;
    double m_timeToExpiry;
    bool m_fixBeta;
};

// Calculate implied volatility using Hagan's formula
double QuantLibSABRModel::impliedVolatility(
    double forward,
    double strike,
    double timeToExpiry) const {

    try {
        // Use QuantLib's SABR implementation
        return QuantLib::sabrVolatility(strike, forward, timeToExpiry, m_alpha, m_beta, m_nu, m_rho);
    } catch (std::exception& e) {
        std::cerr << "Error in SABR implied volatility calculation: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in SABR implied volatility calculation" << std::endl;
        return -1.0;
    }
}

// Calculate implied volatility at a vector of strikes
std::vector<double> QuantLibSABRModel::impliedVolatilityCurve(
    double forward,
    const std::vector<double>& strikes,
    double timeToExpiry) const {

    std::vector<double> vols;
    vols.reserve(strikes.size());

    for (double strike : strikes) {
        vols.push_back(impliedVolatility(forward, strike, timeToExpiry));
    }

    return vols;
}

// Calculate option price using SABR model
double QuantLibSABRModel::optionPrice(
    double forward,
    double strike,
    double timeToExpiry,
    double discountFactor,
    QuantLib::Option::Type optionType) const {

    try {
        // First calculate the implied volatility using SABR model
        double vol = impliedVolatility(forward, strike, timeToExpiry);

        // Then use Black formula to calculate the option price
        QuantLib::BlackCalculator black(
            QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff>(
                new QuantLib::PlainVanillaPayoff(optionType, strike)),
            forward,
            vol * std::sqrt(timeToExpiry),
            discountFactor);

        return black.value();

    } catch (std::exception& e) {
        std::cerr << "Error in SABR option pricing: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error in SABR option pricing" << std::endl;
        return -1.0;
    }
}

// Calculate option Greeks using SABR model
QuantLib::Greeks QuantLibSABRModel::optionGreeks(
    double forward,
    double strike,
    double timeToExpiry,
    double discountFactor,
    QuantLib::Option::Type optionType) const {

    try {
        // First calculate the implied volatility using SABR model
        double vol = impliedVolatility(forward, strike, timeToExpiry);

        // Then use Black formula to calculate the Greeks
        QuantLib::BlackCalculator black(
            QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff>(
                new QuantLib::PlainVanillaPayoff(optionType, strike)),
            forward,
            vol * std::sqrt(timeToExpiry),
            discountFactor);

        // Calculate Greeks
        QuantLib::Greeks greeks;

        // Delta and Gamma with respect to the forward (not the spot)
        greeks.delta = black.delta(forward);
        greeks.gamma = black.gamma(forward);

        // Theta: need to take derivative of price with respect to time
        // Approximate by calculating price at a slightly different time
        const double eps = 1.0/365.0; // One day
        double volEps = impliedVolatility(forward, strike, timeToExpiry - eps);

        QuantLib::BlackCalculator blackEps(
            QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff>(
                new QuantLib::PlainVanillaPayoff(optionType, strike)),
            forward,
            volEps * std::sqrt(timeToExpiry - eps),
            discountFactor);

        greeks.theta = (blackEps.value() - black.value()) / eps;

        // Vega with respect to volatility
        greeks.vega = black.vega(timeToExpiry);

        // Rho with respect to interest rate - approximate as derivative of price with respect to discount factor
        const double epsDf = 0.0001;
        QuantLib::BlackCalculator blackDf(
            QuantLib::ext::shared_ptr<QuantLib::StrikedTypePayoff>(
                new QuantLib::PlainVanillaPayoff(optionType, strike)),
            forward,
            vol * std::sqrt(timeToExpiry),
            discountFactor + epsDf);

        greeks.rho = (blackDf.value() - black.value()) / epsDf;

        return greeks;

    } catch (std::exception& e) {
        std::cerr << "Error in SABR Greeks calculation: " << e.what() << std::endl;
        return QuantLib::Greeks();
    } catch (...) {
        std::cerr << "Unknown error in SABR Greeks calculation" << std::endl;
        return QuantLib::Greeks();
    }
}

// Calibrate SABR parameters to market volatilities
void QuantLibSABRModel::calibrate(
    double forward,
    const std::vector<double>& strikes,
    const std::vector<double>& marketVols,
    double timeToExpiry) {

    try {
        // Check input validity
        if (strikes.size() != marketVols.size() || strikes.empty()) {
            throw std::invalid_argument("Strikes and market volatilities must have same non-zero size");
        }

        bool fixBeta = (strikes.size() < 4); // If we have few points, fix beta

        // Set up the optimization problem
        QuantLib::ext::shared_ptr<QuantLib::CostFunction> costFunction(
            new SABRError(forward, strikes, marketVols, timeToExpiry, fixBeta));

        // Initial guess for the SABR parameters
        QuantLib::Array initialValues(fixBeta ? 3 : 4);
        initialValues[0] = m_alpha;

        if (fixBeta) {
            initialValues[1] = m_nu;
            initialValues[2] = m_rho;
        } else {
            initialValues[1] = m_beta;
            initialValues[2] = m_nu;
            initialValues[3] = m_rho;
        }

        // Set up optimization method (Simplex)
        QuantLib::Simplex optimizationMethod(0.01); // Step size

        // Set up end criteria
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Fix: Create a constraint object first, then pass it by reference to Problem constructor
        QuantLib::NoConstraint noConstraint;
        QuantLib::Problem problem(*costFunction, noConstraint, initialValues);

        optimizationMethod.minimize(problem, endCriteria);

        // Get the optimal parameters
        QuantLib::Array xMin = problem.currentValue();

        // Update SABR parameters
        m_alpha = xMin[0];

        if (fixBeta) {
            // Beta is fixed at 0.5
            m_beta = 0.5;
            m_nu = xMin[1];
            m_rho = xMin[2];
        } else {
            m_beta = xMin[1];
            m_nu = xMin[2];
            m_rho = xMin[3];
        }

    } catch (std::exception& e) {
        std::cerr << "Error in SABR calibration: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error in SABR calibration" << std::endl;
    }
}

} // namespace quant