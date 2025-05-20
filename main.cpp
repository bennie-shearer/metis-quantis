#include <iostream>
#include <thread>
#include <chrono>
#include "metis_calculator_api.hpp"
#include "metis_yield_curve_api.hpp"
#include "quantlib_calculator_api.hpp"
#include "metis_webserver.hpp"
#include "metis_json.hpp"

int main() {
    try {
        // Define server port
        int port = 8080;
        simple_web::SimpleWebServer server(port);

        // Register APIs using functions from metis_calculator_api.cpp
        registerQuantLibCalculatorAPI(server);
        registerYieldCurveAPI(server);

        // Use the namespaced functions if needed
        yield_curve::setupYieldCurveEndpoints(server);
        quantlib_calculator::setupQuantLibEndpoints(server);

        std::cout << "Starting server on port " << port << "..." << std::endl;
        server.start();

        // Display the URL for easy access
        std::cout << "Server is running at: http://localhost:" << port << std::endl;
        std::cout << "API endpoints available at:" << std::endl;
        std::cout << "  - http://localhost:" << port << "/api/calculator" << std::endl;
        std::cout << "  - http://localhost:" << port << "/api/yield-curve" << std::endl;

        // Keep the server running until user presses Enter
        std::cout << "\nPress Enter to stop the server." << std::endl;
        std::cin.get();

        // Stop the server gracefully
        server.stop();

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}