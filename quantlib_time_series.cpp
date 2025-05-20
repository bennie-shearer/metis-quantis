// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_time_series.cpp

#include "quantlib_time_series.hpp"
#include <cmath>

QuantLibTimeSeries::MAModel QuantLibTimeSeries::estimateMA(int q) {
    try {
        // Check that we have enough data
        if (m_values.size() <= static_cast<size_t>(q)) {
            throw std::runtime_error("Need more data points than MA order");
        }

        // MA models are typically estimated using maximum likelihood
        // This is a simplified implementation using method of moments

        // Create result
        MAModel model;
        model.coefficients.resize(q, 0.0);

        // Prepare data
        size_t n = m_values.size();
        double mean = std::accumulate(m_values.begin(), m_values.end(), 0.0) / n;

        // Set constant term
        model.constant = mean;

        // Subtract mean
        std::vector<double> centered(n);
        for (size_t i = 0; i < n; ++i) {
            centered[i] = m_values[i] - mean;
        }

        // Calculate sample autocorrelation
        std::vector<double> autocorr = calculateAutocorrelation(q + 1);

        // For MA(q) model, the theoretical autocorrelation is:
        // rho(k) = (theta_k + theta_1 * theta_{k+1} + ... + theta_{q-k} * theta_q) / (1 + theta_1^2 + ... + theta_q^2)
        // for k = 1, 2, ..., q, and rho(k) = 0 for k > q

        // Set up the optimization problem
        class MACostFunction : public QuantLib::CostFunction {
        private:
            std::vector<double> sampleAutocorr_;
            int q_;

        public:
            MACostFunction(const std::vector<double>& autocorr, int q)
                : sampleAutocorr_(autocorr), q_(q) {}

            QuantLib::Real value(const QuantLib::Array& params) const override {
                // Calculate theoretical autocorrelation
                std::vector<double> theta(q_, 0.0);
                for (int i = 0; i < q_; ++i) {
                    theta[i] = params[i];
                }

                // Calculate denominator
                double denom = 1.0;
                for (int i = 0; i < q_; ++i) {
                    denom += theta[i] * theta[i];
                }

                // Calculate theoretical autocorrelation
                std::vector<double> theorAutocorr(q_ + 1, 0.0);
                for (int k = 1; k <= q_; ++k) {
                    double num = 0.0;
                    for (int j = 0; j < q_ - k; ++j) {
                        num += theta[j] * theta[j + k];
                    }
                    theorAutocorr[k] = num / denom;
                }

                // Calculate sum of squared errors
                double sse = 0.0;
                for (int k = 1; k <= q_; ++k) {
                    double diff = sampleAutocorr_[k] - theorAutocorr[k];
                    sse += diff * diff;
                }

                return sse;
            }

            // Add the missing method
            QuantLib::Array values(const QuantLib::Array& params) const override {
                QuantLib::Array result(1);
                result[0] = value(params);
                return result;
            }
        };

        // Create initial parameters
        QuantLib::Array initial(q, 0.1);

        // Create a NoConstraint object first
        QuantLib::NoConstraint noConstraint;

        // Set up the optimization problem with proper constraint handling
        MACostFunction costFunction(autocorr, q);
        QuantLib::Problem problem(costFunction, noConstraint, initial);

        // Set up the optimizer
        QuantLib::Simplex optimizer(0.1);

        // Set up the end criteria
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Optimize
        optimizer.minimize(problem, endCriteria);

        // Get the optimized parameters
        QuantLib::Array theta = problem.currentValue();

        // Fill the model coefficients
        for (int i = 0; i < q; ++i) {
            model.coefficients[i] = theta[i];
        }

        // Calculate innovation variance (simplified)
        double gamma0 = 0.0;
        for (double x : centered) {
            gamma0 += x * x;
        }
        gamma0 /= n;

        double denom = 1.0;
        for (int i = 0; i < q; ++i) {
            denom += theta[i] * theta[i];
        }

        model.sigma = std::sqrt(gamma0 / denom);

        return model;

    } catch (std::exception& e) {
        std::cerr << "Error estimating MA model: " << e.what() << std::endl;
        return MAModel();
    } catch (...) {
        std::cerr << "Unknown error estimating MA model" << std::endl;
        return MAModel();
    }
}

QuantLibTimeSeries::GARCHModel QuantLibTimeSeries::estimateGARCH(int p, int q) {
    try {
        // Check that we have enough data
        if (m_values.size() <= static_cast<size_t>(p + q + 1)) {
            throw std::runtime_error("Need more data points than GARCH order");
        }

        // Simplified GARCH(1,1) implementation
        if (p != 1 || q != 1) {
            throw std::runtime_error("Only GARCH(1,1) is currently supported");
        }

        // Create result
        GARCHModel model;

        // Prepare data
        size_t n = m_values.size();
        double mean = std::accumulate(m_values.begin(), m_values.end(), 0.0) / n;

        // Calculate residuals
        std::vector<double> residuals(n);
        for (size_t i = 0; i < n; ++i) {
            residuals[i] = m_values[i] - mean;
        }

        // Initialize parameters
        model.omega = 0.01;
        model.alpha = 0.1;
        model.beta = 0.8;

        // Initial standard deviation
        double variance = 0.0;
        for (double r : residuals) {
            variance += r * r;
        }
        variance /= n;
        model.sigma = std::sqrt(variance);

        // Set up the optimization problem
        class GARCHCostFunction : public QuantLib::CostFunction {
        private:
            std::vector<double> residuals_;

        public:
            GARCHCostFunction(const std::vector<double>& residuals)
                : residuals_(residuals) {}

            QuantLib::Real value(const QuantLib::Array& params) const override {
                // Extract parameters
                double omega = params[0];
                double alpha = params[1];
                double beta = params[2];

                // Check parameter constraints
                if (omega <= 0.0 || alpha <= 0.0 || beta <= 0.0 || alpha + beta >= 1.0) {
                    return 1.0e10; // Large penalty for invalid parameters
                }

                // Calculate log-likelihood
                double logL = 0.0;
                double h = omega / (1.0 - alpha - beta); // Unconditional variance

                for (size_t t = 1; t < residuals_.size(); ++t) {
                    h = omega + alpha * residuals_[t-1] * residuals_[t-1] + beta * h;
                    logL += -0.5 * (std::log(2.0 * M_PI) + std::log(h) + residuals_[t] * residuals_[t] / h);
                }

                return -logL; // Minimize negative log-likelihood
            }

            // Add the missing method
            QuantLib::Array values(const QuantLib::Array& params) const override {
                QuantLib::Array result(1);
                result[0] = value(params);
                return result;
            }
        };

        // Create initial parameters
        QuantLib::Array initial(3);
        initial[0] = model.omega;
        initial[1] = model.alpha;
        initial[2] = model.beta;

        // Create a NoConstraint object first
        QuantLib::NoConstraint noConstraint;

        // Set up the optimization problem with proper constraint handling
        GARCHCostFunction costFunction(residuals);
        QuantLib::Problem problem(costFunction, noConstraint, initial);

        // Set up the optimizer
        QuantLib::LevenbergMarquardt optimizer;

        // Set up the end criteria
        QuantLib::EndCriteria endCriteria(10000, 100, 1e-8, 1e-8, 1e-8);

        // Optimize
        optimizer.minimize(problem, endCriteria);

        // Get the optimized parameters
        QuantLib::Array optimal = problem.currentValue();

        // Fill the model parameters
        model.omega = optimal[0];
        model.alpha = optimal[1];
        model.beta = optimal[2];

        return model;

    } catch (std::exception& e) {
        std::cerr << "Error estimating GARCH model: " << e.what() << std::endl;
        return GARCHModel();
    } catch (...) {
        std::cerr << "Unknown error estimating GARCH model" << std::endl;
        return GARCHModel();
    }
}