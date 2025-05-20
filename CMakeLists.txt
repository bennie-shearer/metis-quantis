cmake_minimum_required(VERSION 3.10)
project(quant-boost)

set(CMAKE_CXX_STANDARD 20)
message(STATUS "Using C++20 standard")

# Find Boost
set(BOOST_ROOT "C:/Boost/boost_1_82_0")
find_package(Boost 1.82.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")

# QuantLib setup
set(QUANTLIB_INCLUDE_DIR "C:/quantlib/QuantLib-1.30/QuantLib-1.30")
include_directories(${QUANTLIB_INCLUDE_DIR})
message(STATUS "QuantLib include dir: ${QUANTLIB_INCLUDE_DIR}")

# Set the path to the QuantLib static library
set(QUANTLIB_LIBRARY "${CMAKE_BINARY_DIR}/libQuantLib.a")
message(STATUS "Using QuantLib library: ${QUANTLIB_LIBRARY}")

# Add socket library for Windows
if(WIN32)
    set(SOCKET_LIBRARY ws2_32)
endif()

# Define all library targets
add_library(metis_webserver_lib STATIC metis_webserver.cpp)
add_library(metis_calculator_api_lib STATIC metis_calculator_api.cpp)
add_library(metis_yield_curve_api_lib STATIC metis_yield_curve_api.cpp)
add_library(quantlib_bond_calculator_lib STATIC quantlib_bond_calculator.cpp)
add_library(quantlib_binomial_tree_lib STATIC quantlib_binomial_tree.cpp)
add_library(quantlib_bonds_lib STATIC quantlib_bonds.cpp)
add_library(quantlib_caps_floors_lib STATIC quantlib_caps_floors.cpp)
add_library(quantlib_calculator_api_lib STATIC quantlib_calculator_api.cpp)
add_library(quantlib_commodity_lib STATIC quantlib_commodity.cpp)
add_library(quantlib_caps_floors_api_lib STATIC quantlib_caps_floors_api.cpp)

# Make sure all libraries have access to the project include directory
target_include_directories(metis_calculator_api_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(metis_yield_curve_api_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_calculator_api_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_bond_calculator_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_bonds_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_caps_floors_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_caps_floors_api_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_binomial_tree_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantlib_commodity_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Define the main executable
add_executable(main main.cpp)

# Add linker flags to fix duplicate sections issue
if(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
endif()

# Link the main executable with all libraries
target_link_libraries(main
        metis_webserver_lib
        metis_calculator_api_lib
        metis_yield_curve_api_lib
        quantlib_bond_calculator_lib
        quantlib_binomial_tree_lib
        quantlib_bonds_lib
        quantlib_caps_floors_lib
        quantlib_calculator_api_lib
        quantlib_commodity_lib
        quantlib_caps_floors_api_lib
        ${QUANTLIB_LIBRARY}
        ${SOCKET_LIBRARY}
)