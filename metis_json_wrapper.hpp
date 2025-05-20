// File: C:/Users/bshearer/CLionProjects/quant-boost/metis_json_wrapper.hpp
#ifndef METIS_JSON_WRAPPER_H
#define METIS_JSON_WRAPPER_H

#include "metis_json.hpp"

namespace quant {

// Helper wrapper class with const-correct methods for SimpleJSON
class JsonWrapper {
public:
    JsonWrapper(const simple_json::SimpleJSON& json) : m_json(json) {}

    // Access methods for different types
    std::string asString() const {
        // Create a non-const copy for access
        simple_json::SimpleJSON copy = m_json;
        return copy.asString();
    }

    double asDouble() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.asDouble();
    }

    int asInt() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.asInt();
    }

    bool asBool() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.asBool();
    }

    std::vector<simple_json::SimpleJSON> asArray() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.asArray();
    }

    // Check type methods
    bool isString() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isString();
    }

    bool isNumber() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isNumber();
    }

    bool isInt() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isInt();
    }

    bool isBool() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isBool();
    }

    bool isArray() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isArray();
    }

    bool isObject() const {
        simple_json::SimpleJSON copy = m_json;
        return copy.isObject();
    }

    // Access an element by key (for objects)
    JsonWrapper operator[](const std::string& key) const {
        simple_json::SimpleJSON copy = m_json;
        return JsonWrapper(copy[key]);
    }

    // Access an element by index (for arrays)
    JsonWrapper operator[](size_t index) const {
        simple_json::SimpleJSON copy = m_json;
        return JsonWrapper(copy[index]);
    }

private:
    simple_json::SimpleJSON m_json;
};

} // namespace quant

#endif // METIS_JSON_WRAPPER_H