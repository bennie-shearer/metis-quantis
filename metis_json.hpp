#ifndef METIS_JSON_HPP
#define METIS_JSON_HPP

#include <string>
#include <map>
#include <vector>
#include <variant>

namespace simple_json {
    class SimpleJSON {
    private:
        // Internal storage for various value types
        std::map<std::string, std::variant<std::string, double, bool, int>> values;

    public:
        // Constructors
        SimpleJSON() = default;

        // Methods to set values
        void setString(const std::string& key, const std::string& value) {
            values[key] = value;
        }

        void setNumber(const std::string& key, double value) {
            values[key] = value;
        }

        void setBool(const std::string& key, bool value) {
            values[key] = value;
        }

        // Methods to get values
        std::string getString(const std::string& key) const {
            try {
                return std::get<std::string>(values.at(key));
            } catch (...) {
                return "";
            }
        }

        double getNumber(const std::string& key) const {
            try {
                return std::get<double>(values.at(key));
            } catch (...) {
                return 0.0;
            }
        }

        bool getBool(const std::string& key) const {
            try {
                return std::get<bool>(values.at(key));
            } catch (...) {
                return false;
            }
        }

        // Method to convert to string
        std::string toString() const {
            std::string result = "{";
            bool first = true;

            for (const auto& [key, value] : values) {
                if (!first) result += ",";
                first = false;

                result += "\"" + key + "\":";

                if (std::holds_alternative<std::string>(value)) {
                    result += "\"" + std::get<std::string>(value) + "\"";
                } else if (std::holds_alternative<double>(value)) {
                    result += std::to_string(std::get<double>(value));
                } else if (std::holds_alternative<bool>(value)) {
                    result += std::get<bool>(value) ? "true" : "false";
                } else if (std::holds_alternative<int>(value)) {
                    result += std::to_string(std::get<int>(value));
                }
            }

            result += "}";
            return result;
        }

        // Add array support (simplified)
        void setArray(const std::string& key, const SimpleJSON& arrayJson) {
            // For simplicity, we just store the JSON string representation
            setString(key, arrayJson.toString());
        }
    };
}

#endif // METIS_JSON_HPP