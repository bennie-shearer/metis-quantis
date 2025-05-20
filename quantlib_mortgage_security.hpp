// File: C:/Users/bshearer/CLionProjects/quant-boost/quantlib_mortgage_security.hpp
#ifndef QUANTLIB_MORTGAGE_SECURITY_HPP
#define QUANTLIB_MORTGAGE_SECURITY_HPP

#include <vector>
#include <memory>
#include <ql/quantlib.hpp>

namespace quant {

class QuantLibMortgageSecurity {
public:
    // Struct for mortgage pool parameters
    struct MortgagePoolParams {
        double notional;
        double wac;          // Weighted Average Coupon
        double wam;          // Weighted Average Maturity (in years)
        double psa;          // PSA prepayment assumption (100 = 100% PSA)
        QuantLib::DayCounter dayCounter;
    };

    // Struct for mortgage tranche parameters
    struct TrancheParams {
        double principal;
        double coupon;
        int priority;        // Priority in cash flow waterfall (lower numbers get paid first)
        bool isInterestPaying;
        bool isPrincipalPaying;
    };

    // Constructor
    QuantLibMortgageSecurity(
        const MortgagePoolParams& poolParams = {1'000'000, 0.045, 30.0, 100.0, QuantLib::Thirty360(QuantLib::Thirty360::BondBasis)});

    // Add a tranche to the mortgage security
    void addTranche(const TrancheParams& trancheParams);

    // Calculate pool cash flows
    std::pair<std::vector<double>, std::vector<double>> calculatePoolCashFlows(
        const QuantLib::Date& startDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Calculate tranche cash flows
    std::vector<std::pair<std::vector<double>, std::vector<double>>> calculateTrancheCashFlows(
        const QuantLib::Date& startDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Calculate tranche value
    double calculateTrancheValue(
        size_t trancheIndex,
        const QuantLib::Date& valuationDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve);

    // Calculate tranche yield
    double calculateTrancheYield(
        size_t trancheIndex,
        double price,
        const QuantLib::Date& startDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

    // Calculate tranche weighted average life
    double calculateTrancheWAL(
        size_t trancheIndex,
        const QuantLib::Date& startDate,
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure);

private:
    // Calculate prepayment rate based on PSA model
    double calculatePrepaymentRate(double age, double psa);

    // Pool parameters
    MortgagePoolParams m_poolParams;

    // Tranches
    std::vector<TrancheParams> m_tranches;
};

} // namespace quant

#endif // QUANTLIB_MORTGAGE_SECURITY_HPP