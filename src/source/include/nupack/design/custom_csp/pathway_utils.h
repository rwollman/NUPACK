#pragma once

#include "pair_probabilities.h"

#include "types.h"

#include "adapter.h"

#include <vector>
#include <sstream>
#include <type_traits>

namespace nupack {
namespace custom_csp {

int sample_weighted_int(const vec<real> & weights);
real get_current_time();

/*
 * This format is used for partition functions and concentrations,
 * anything that needs to be a scale-free parameter
 */
struct scale_format_t {
    friend std::ostream & operator<<(std::ostream & out, scale_format_t const &) {
        out << std::setw(6) << std::setprecision(3) << std::fixed;
        return out;
    }
};

struct exp_format_t {
    friend std::ostream & operator<<(std::ostream & out, exp_format_t const &) {
        out << std::setw(10) << std::setprecision(6) << std::scientific;
        return out;
    }
};

struct flt_format_t {
    friend std::ostream & operator<<(std::ostream & out, flt_format_t const &) {
        out << std::setw(10); out << std::setprecision(6); out << std::fixed;
        return out;
    }
};

struct longflt_format_t {
    friend std::ostream & operator<<(std::ostream & out, longflt_format_t const &) {
        out << std::setw(14) << std::setprecision(10) << std::fixed;
        return out;
    }
};

struct bool_format_t {
    friend std::ostream & operator<<(std::ostream & out, bool_format_t const &) {
        out << std::boolalpha;
        return out;
    }
};

struct param_format_t {
    friend std::ostream & operator<<(std::ostream & out, param_format_t const &) {
        out << std::setw(7);
        return out;
    }
};

struct null_format_t {
    friend std::ostream & operator<<(std::ostream & out, null_format_t const &) {return out;}
};

static constexpr scale_format_t scale_format = {};
static constexpr exp_format_t exp_format = {};
static constexpr flt_format_t flt_format = {};
static constexpr longflt_format_t longflt_format = {};
static constexpr bool_format_t bool_format = {};
static constexpr param_format_t param_format = {};
static constexpr null_format_t null_format = {};

template < class T, int_if < !is_iterable<T> && !is_pair<T> > = 0 >
string to_string(T t);

template <class T, int_if<is_pair<T>> = 0>
string to_string(T t);

template <class T, int_if<is_iterable<T>> = 0>
string to_string(T t);

// string to_string(bool b);
string to_string(string b);

template < class T, int_if < !is_iterable<T> && !is_pair<T> >>
inline string to_string(T t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <class T, int_if<is_pair<T>>>
inline string to_string(T t) {
    std::stringstream ss;
    ss << '{' << to_string(t.first) << ": " << to_string(t.second) << "}";
    return ss.str();
}

template <class T, int_if<is_iterable<T>>>
inline string to_string(T t) {
    std::stringstream ss;
    ss << '(';
    auto i = t.size();
    for (auto & el : t) {
        ss << (--i == 0 ? to_string(el) : to_string(el) + ", ");
    }
    ss << ')';
    return ss.str();
}

inline string to_string(bool b) {
    return b ? "true" : "false";
};

inline string to_string(string str) {
    return str;
};

}
}

