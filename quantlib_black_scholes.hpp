// File: metis-quantis/ql/pricingengines/blackcommodityengine.hpp

#ifndef quantlib_black_commodity_engine_hpp
#define quantlib_black_commodity_engine_hpp

#include <ql/instruments/commodity.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    class BlackCommodityOptionEngine : public CommodityOption::engine {
      public:
        BlackCommodityOptionEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
        void calculate() const override;

      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
    };

}

#endif

// File: metis-quantis/ql/pricingengines/blackcommodityengine.cpp

#include <ql/pricingengines/blackcommodityengine.hpp>
#include <ql/pricingengines/blackcalculator.hpp>
#include <ql/exercise.hpp>

namespace QuantLib {

    BlackCommodityOptionEngine::BlackCommodityOptionEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process)
    : process_(process) {
        registerWith(process_);
    }

    void BlackCommodityOptionEngine::calculate() const {
        // Extract option data
        auto payoff = ext::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-striked payoff given");

        Date exerciseDate = arguments_.exercise->lastDate();
        auto forward = arguments_.forward;

        // Get relevant market data
        Real spotPrice = forward->spotPrice();
        Date settlementDate = process_->riskFreeRate()->referenceDate();
        Time t = process_->time(exerciseDate);

        // Early exercise handling
        if (arguments_.exercise->type() == Exercise::European) {
            // Collect inputs for Black formula
            Real strike = payoff->strike();
            Real forwardPrice = spotPrice * process_->dividendYield()->discount(t) /
                                process_->riskFreeRate()->discount(t);
            Real stdDev = process_->blackVolatility()->blackVol(exerciseDate, strike) *
                         std::sqrt(t);
            DiscountFactor discount = process_->riskFreeRate()->discount(t);

            // Option price calculation using Black formula
            BlackCalculator black(payoff, forwardPrice, stdDev, discount);
            results_.value = black.value();
            results_.delta = black.delta(spotPrice);
            results_.gamma = black.gamma(spotPrice);
            results_.theta = black.theta(spotPrice, t);
            results_.vega = black.vega(t);
            results_.rho = black.rho(t);
        } else {
            QL_FAIL("non-European exercise not supported");
        }
    }
}