// pyfmt.hpp - Python-compatible float formatting.
//
// The Python code builds output directory names out of f-strings such as
//   f"delta-max_{args.delta_max}"      -> str(float)  == repr(float)
//   f"delta_min_{args.tolerance:.1e}"  -> '1.0e-10'
//   f"dt_{env.dt:.3f}"                 -> '0.300'
// Reproducing those byte-for-byte keeps the C++/Rust runs writing into the same
// results directories as the Python reference implementation.
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace dhj {

// Equivalent of Python's repr(float) / str(float): shortest representation that
// round-trips, switching to scientific notation when exp < -4 or exp >= 16.
inline std::string py_repr(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    if (v == 0.0) return std::signbit(v) ? "-0.0" : "0.0";

    char buf[64];
    int prec = 0;
    for (; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof buf, "%.*e", prec, v);
        if (std::strtod(buf, nullptr) == v) break;
    }

    std::string s(buf);
    const std::size_t epos = s.find('e');
    std::string mant = s.substr(0, epos);
    const int exp = std::atoi(s.c_str() + epos + 1);

    bool neg = false;
    if (!mant.empty() && mant[0] == '-') { neg = true; mant.erase(0, 1); }

    std::string digits;
    for (char c : mant) if (c != '.') digits += c;
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    const int ndig = static_cast<int>(digits.size());

    std::string out;
    if (exp < -4 || exp >= 16) {
        out = digits.substr(0, 1);
        if (ndig > 1) out += "." + digits.substr(1);
        char eb[16];
        std::snprintf(eb, sizeof eb, "e%c%02d", exp < 0 ? '-' : '+', std::abs(exp));
        out += eb;
    } else if (exp >= 0) {
        if (ndig > exp + 1) out = digits.substr(0, exp + 1) + "." + digits.substr(exp + 1);
        else                out = digits + std::string(exp + 1 - ndig, '0') + ".0";
    } else {
        out = "0." + std::string(-exp - 1, '0') + digits;
    }
    return neg ? "-" + out : out;
}

// Python's format(v, '.<prec>e'); C's %e already matches (>= 2 exponent digits).
inline std::string py_exp(double v, int prec) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", prec, v);
    return std::string(buf);
}

// Python's format(v, '.<prec>f').
inline std::string py_fixed(double v, int prec) {
    char buf[512];
    std::snprintf(buf, sizeof buf, "%.*f", prec, v);
    return std::string(buf);
}

inline std::string py_bool(bool b) { return b ? "True" : "False"; }

}  // namespace dhj
