// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_mortgage_security.cpp

#include "quantlib_mortgage_security.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace quant {

// Constructor implementation
QuantLibMortgageSecurity::QuantLibMortgageSecurity(const MortgagePoolParams& poolParams)
    : m_poolParams(poolParams) {
}

void QuantLibMortgageSecurity::addTranche(const TrancheParams& trancheParams) {
    m_tranches.push_back(trancheParams);
}

std::pair<std::vector<double>, std::vector<double>> QuantLibMortgageSecurity::calculatePoolCashFlows(
    const QuantLib::Date& startDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) {

    try {
        // Calculate total number of months
        int totalMonths = static_cast<int>(m_poolParams.wam * 12);

        // Initialize remaining balance, principal and interest cash flow vectors
        double remainingBalance = m_poolParams.notional;
        std::vector<double> principalCashFlows(totalMonths, 0.0);
        std::vector<double> interestCashFlows(totalMonths, 0.0);

        // Calculate monthly payment (mortgage formula: P = L[c(1+c)^n]/[(1+c)^n-1])
        double monthlyRate = m_poolParams.wac / 12.0;
        double payment = m_poolParams.notional * monthlyRate *
                         std::pow(1.0 + monthlyRate, totalMonths) /
                         (std::pow(1.0 + monthlyRate, totalMonths) - 1.0);

        // For each month, calculate scheduled principal, interest, and prepayment
        for (int month = 0; month < totalMonths; ++month) {
            // Skip if no balance remains
            if (remainingBalance <= 0.0) break;

            // Calculate scheduled interest
            double interest = remainingBalance * monthlyRate;
            interestCashFlows[month] += interest;

            // Calculate scheduled principal
            double scheduledPrincipal = payment - interest;

            // Ensure we don't pay more than the remaining balance
            scheduledPrincipal = std::min(scheduledPrincipal, remainingBalance);

            // Calculate prepayment (CPR based on PSA)
            double cpr = calculatePrepaymentRate(month + 1, m_poolParams.psa);
            double smm = 1.0 - std::pow(1.0 - cpr, 1.0/12.0); // Convert annual CPR to monthly SMM
            double prepayment = (remainingBalance - scheduledPrincipal) * smm;

            // Total principal payment
            double totalPrincipal = scheduledPrincipal + prepayment;
            principalCashFlows[month] += totalPrincipal;

            // Update remaining balance
            remainingBalance -= totalPrincipal;
        }

        return {principalCashFlows, interestCashFlows};
    } catch (std::exception& e) {
        std::cerr << "Error calculating pool cash flows: " << e.what() << std::endl;
        return {std::vector<double>(), std::vector<double>()};
    } catch (...) {
        std::cerr << "Unknown error calculating pool cash flows" << std::endl;
        return {std::vector<double>(), std::vector<double>()};
    }
}

double QuantLibMortgageSecurity::calculateTrancheValue(
    size_t trancheIndex,
    const QuantLib::Date& valuationDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve) {

    try {
        // Verify tranche index
        if (trancheIndex >= m_tranches.size()) {
            throw std::out_of_range("Tranche index out of range");
        }

        // Get tranche cash flows
        auto trancheCashFlows = calculateTrancheCashFlows(valuationDate, discountCurve);
        auto [principalCFs, interestCFs] = trancheCashFlows[trancheIndex];

        // Calculate number of months
        int numMonths = static_cast<int>(principalCFs.size());

        // Calculate present value of all cash flows
        double pv = 0.0;
        for (int month = 0; month < numMonths; ++month) {
            // Calculate payment date
            QuantLib::Date paymentDate = valuationDate + QuantLib::Period(month + 1, QuantLib::Months);

            // Calculate discount factor
            double df = discountCurve->discount(paymentDate);

            // Add discounted cash flow to present value
            pv += (principalCFs[month] + interestCFs[month]) * df;
        }

        return pv;
    } catch (std::exception& e) {
        std::cerr << "Error calculating tranche value: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating tranche value" << std::endl;
        return -1.0;
    }
}

double QuantLibMortgageSecurity::calculateTrancheWAL(
    size_t trancheIndex,
    const QuantLib::Date& startDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) {

    try {
        // Verify tranche index
        if (trancheIndex >= m_tranches.size()) {
            throw std::out_of_range("Tranche index out of range");
        }

        // Get tranche cash flows
        auto trancheCashFlows = calculateTrancheCashFlows(startDate, termStructure);
        auto principalCFs = trancheCashFlows[trancheIndex].first;

        // Calculate total principal and weighted sum
        double totalPrincipal = 0.0;
        double weightedSum = 0.0;

        for (size_t i = 0; i < principalCFs.size(); ++i) {
            double timeInYears = (i + 1) / 12.0; // Convert months to years
            weightedSum += principalCFs[i] * timeInYears;
            totalPrincipal += principalCFs[i];
        }

        // Calculate WAL
        if (totalPrincipal > 0.0) {
            return weightedSum / totalPrincipal;
        } else {
            return 0.0;
        }
    } catch (std::exception& e) {
        std::cerr << "Error calculating tranche WAL: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating tranche WAL" << std::endl;
        return -1.0;
    }
}

double QuantLibMortgageSecurity::calculateTrancheYield(
    size_t trancheIndex,
    double price,
    const QuantLib::Date& startDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) {

    try {
        // Verify tranche index
        if (trancheIndex >= m_tranches.size()) {
            throw std::out_of_range("Tranche index out of range");
        }

        // Get tranche cash flows
        auto trancheCashFlows = calculateTrancheCashFlows(startDate, termStructure);
        auto [principalCFs, interestCFs] = trancheCashFlows[trancheIndex];

        // Combine principal and interest cash flows
        std::vector<double> totalCFs(principalCFs.size());
        for (size_t i = 0; i < principalCFs.size(); ++i) {
            totalCFs[i] = principalCFs[i] + interestCFs[i];
        }

        // Create a set of dates for the cashflows
        std::vector<QuantLib::Date> dates;
        for (size_t i = 0; i < totalCFs.size(); ++i) {
            dates.push_back(startDate + QuantLib::Period(static_cast<int>(i + 1), QuantLib::Months));
        }

        // Create a QuantLib::Leg with the cashflows - using boost::shared_ptr instead of std::shared_ptr
        QuantLib::Leg leg;
        for (size_t i = 0; i < totalCFs.size(); ++i) {
            if (totalCFs[i] > 0.0) {
                leg.push_back(boost::make_shared<QuantLib::SimpleCashFlow>(totalCFs[i], dates[i]));
            }
        }

        // Add the negative initial investment to the cash flows
        leg.insert(leg.begin(), boost::make_shared<QuantLib::SimpleCashFlow>(-price, startDate));

        // Initialize solver parameters
        QuantLib::Real guess = 0.05; // Initial guess for yield
        QuantLib::Size maxIterations = 1000;
        QuantLib::Real accuracy = 1.0e-10;

        // Calculate IRR using QuantLib's internal rate function
        QuantLib::Real result = QuantLib::CashFlows::yield(
            leg,
            price,
            m_poolParams.dayCounter,
            QuantLib::Compounded,
            QuantLib::Monthly,
            false,
            startDate,
            startDate,
            accuracy,
            maxIterations,
            guess);

        // Convert monthly yield to annual yield
        return std::pow(1.0 + result, 12.0) - 1.0;
    } catch (std::exception& e) {
        std::cerr << "Error calculating tranche yield: " << e.what() << std::endl;
        return -1.0;
    } catch (...) {
        std::cerr << "Unknown error calculating tranche yield" << std::endl;
        return -1.0;
    }
}

std::vector<std::pair<std::vector<double>, std::vector<double>>> QuantLibMortgageSecurity::calculateTrancheCashFlows(
    const QuantLib::Date& startDate,
    const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure) {

    try {
        // Calculate pool cash flows
        auto [poolPrincipalCFs, poolInterestCFs] = calculatePoolCashFlows(startDate, termStructure);

        // Initialize tranche cash flows
        std::vector<std::pair<std::vector<double>, std::vector<double>>> trancheCashFlows;
        for (size_t i = 0; i < m_tranches.size(); ++i) {
            trancheCashFlows.push_back({
                std::vector<double>(poolPrincipalCFs.size(), 0.0), // Principal
                std::vector<double>(poolInterestCFs.size(), 0.0)   // Interest
            });
        }

        // Sort tranches by priority
        std::vector<size_t> trancheOrder(m_tranches.size());
        for (size_t i = 0; i < m_tranches.size(); ++i) {
            trancheOrder[i] = i;
        }
        std::sort(trancheOrder.begin(), trancheOrder.end(), [this](size_t a, size_t b) {
            return m_tranches[a].priority < m_tranches[b].priority;
        });

        // Initialize remaining tranche principal
        std::vector<double> tranchePrincipalRemaining(m_tranches.size());
        for (size_t i = 0; i < m_tranches.size(); ++i) {
            tranchePrincipalRemaining[i] = m_tranches[i].principal;
        }

        // Distribute cash flows to tranches for each month
        for (size_t month = 0; month < poolPrincipalCFs.size(); ++month) {
            // Distribute interest
            double interestRemaining = poolInterestCFs[month];
            for (size_t i : trancheOrder) {
                if (m_tranches[i].isInterestPaying && tranchePrincipalRemaining[i] > 0.0) {
                    // Calculate interest based on tranche coupon
                    double trancheInterest = tranchePrincipalRemaining[i] * (m_tranches[i].coupon / 12.0);

                    // Ensure we don't distribute more than available
                    trancheInterest = std::min(trancheInterest, interestRemaining);
                    trancheCashFlows[i].second[month] = trancheInterest;
                    interestRemaining -= trancheInterest;

                    if (interestRemaining <= 0.0) break;
                }
            }

            // Distribute principal
            double principalRemaining = poolPrincipalCFs[month];
            for (size_t i : trancheOrder) {
                if (m_tranches[i].isPrincipalPaying && tranchePrincipalRemaining[i] > 0.0) {
                    // Calculate principal distribution
                    double tranchePrincipal = std::min(tranchePrincipalRemaining[i], principalRemaining);
                    trancheCashFlows[i].first[month] = tranchePrincipal;
                    tranchePrincipalRemaining[i] -= tranchePrincipal;
                    principalRemaining -= tranchePrincipal;

                    if (principalRemaining <= 0.0) break;
                }
            }
        }

        return trancheCashFlows;
    } catch (std::exception& e) {
        std::cerr << "Error calculating tranche cash flows: " << e.what() << std::endl;
        return {};
    } catch (...) {
        std::cerr << "Unknown error calculating tranche cash flows" << std::endl;
        return {};
    }
}

double QuantLibMortgageSecurity::calculatePrepaymentRate(double age, double psa) {
    try {
        // Standard PSA model:
        // CPR = min(age/30, 1) * (PSA/100) * 0.06
        // where age is in months, and 0.06 is the terminal CPR for 100% PSA (6%)

        double standardCPR = 0.06; // 6% annual prepayment rate for 100% PSA
        double rampFactor = std::min(age / 30.0, 1.0); // Linear ramp-up over first 30 months

        return rampFactor * (psa / 100.0) * standardCPR;
    } catch (std::exception& e) {
        std::cerr << "Error calculating prepayment rate: " << e.what() << std::endl;
        return 0.0;
    } catch (...) {
        std::cerr << "Unknown error calculating prepayment rate" << std::endl;
        return 0.0;
    }
}

} // namespace quant