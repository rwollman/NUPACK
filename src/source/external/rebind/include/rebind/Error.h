#pragma once
#include "Signature.h"

#include <stdexcept>
#include <string_view>
#include <vector>
#include <typeindex>
#include <algorithm>
#include <string>
#include <any>
#include <deque>
#include <optional>

namespace rebind {

/******************************************************************************/

/// Exception for something the API user caused
struct ClientError : std::exception {
    std::string_view message;
    explicit ClientError(std::string_view const &s) noexcept : message(s) {}
    char const * what() const noexcept override {return message.empty() ? "rebind::ClientError" : message.data();}
};

/******************************************************************************/

struct DispatchError : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

/******************************************************************************/

/// Exception for wrong number of arguments
struct WrongNumber : DispatchError {
    unsigned int expected, received;
    WrongNumber(unsigned int n0, unsigned int n)
        : DispatchError("wrong number of arguments"), expected(n0), received(n) {}
};

/******************************************************************************/

/// Exception for wrong type of an argument
struct WrongType : DispatchError {
    std::vector<unsigned int> indices;
    std::string source;
    TypeIndex dest;
    int index, expected, received;

    WrongType(std::string const &n, std::vector<unsigned int> &&v,
              std::string &&s, TypeIndex &&d, int i, int e=0, int r=0) noexcept
        : DispatchError(n), indices(std::move(v)), source(std::move(s)),
          dest(std::move(d)), index(i), expected(e), received(r) {}
};

/******************************************************************************/

struct Dispatch {
    std::string scope;
    Caller caller;
    std::deque<std::any> storage; // deque is used so references don't go bad when doing emplace_back()
    std::vector<unsigned int> indices;
    std::string source;
    TypeIndex dest;
    int index = -1, expected = -1, received = -1;

    std::nullopt_t error() noexcept {return std::nullopt;}

    /// Set error information and return std::nullopt for convenience
    std::nullopt_t error(std::string msg) noexcept {
        scope = std::move(msg);
        return std::nullopt;
    }

    /// Set error information and return std::nullopt for convenience
    std::nullopt_t error(TypeIndex d) noexcept {
        dest = std::move(d);
        return std::nullopt;
    }

    /// Set error information and return std::nullopt for convenience
    std::nullopt_t error(std::string msg, TypeIndex d, int e=-1, int r=-1) noexcept {
        scope = std::move(msg);
        dest = std::move(d);
        expected = e;
        received = r;
        return std::nullopt;
    }

    /// Create exception from the current scopes and messages
    WrongType exception() && noexcept {
        return {std::move(scope), std::move(indices), std::move(source), std::move(dest), index, expected, received};
    }

    /// Store a value which will last the lifetime of a conversion request. Return its address
    template <class T>
    unqualified<T> * store(T &&t) {
        return std::addressof(storage.emplace_back().emplace<unqualified<T>>(static_cast<T &&>(t)));
    }

    Dispatch(Caller c={}, char const *s="mismatched type") : scope(s), caller(std::move(c)) {indices.reserve(8);}
};

/******************************************************************************/

}
