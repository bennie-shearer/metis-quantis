// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_time_series.hpp
#ifndef QUANTLIB_TIME_SERIES_HPP
#define QUANTLIB_TIME_SERIES_HPP

#include <vector>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/constraint.hpp>

class QuantLibTimeSeries {
public:
    // Result structures
    struct MAModel {
        double constant = 0.0;
        std::vector<double> coefficients;
        double sigma = 0.0;
    };

    struct GARCHModel {
        double omega = 0.01;
        double alpha = 0.1;
        double beta = 0.8;
        double sigma = 0.0;
    };

    // Constructor
    QuantLibTimeSeries(const std::vector<double>& values) : m_values(values) {}

    // Model estimation methods
    MAModel estimateMA(int q);
    GARCHModel estimateGARCH(int p, int q);

    // Helper function for autocorrelation calculation
    std::vector<double> calculateAutocorrelation(size_t maxLag) {
        size_t n = m_values.size();
        double mean = std::accumulate(m_values.begin(), m_values.end(), 0.0) / n;

        std::vector<double> centered(n);
        for (size_t i = 0; i < n; ++i) {
            centered[i] = m_values[i] - mean;
        }

        double variance = 0.0;
        for (double c : centered) {
            variance += c * c;
        }
        variance /= n;

        std::vector<double> autocorr(maxLag + 1, 0.0);
        autocorr[0] = 1.0;  // Autocorrelation at lag 0 is 1

        for (size_t lag = 1; lag <= maxLag; ++lag) {
            double sum = 0.0;
            for (size_t i = lag; i < n; ++i) {
                sum += centered[i] * centered[i - lag];
            }
            autocorr[lag] = sum / (n - lag) / variance;
        }

        return autocorr;
    }

private:
    std::vector<double> m_values;
};

#endif // QUANTLIB_TIME_SERIES_HPP