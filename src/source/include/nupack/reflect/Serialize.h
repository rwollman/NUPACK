#pragma once
#include <nlohmann/json.hpp>

#include "../standard/Variant.h"
#include "../standard/Map.h"
#include "../standard/Ptr.h"
#include "../standard/Optional.h"
#include "../reflect/Reflection.h"

/******************************************************************************************/

namespace nupack {

using json = nlohmann::json;

inline json merge_json(json out, json updates) {
    out.insert(updates.begin(), updates.end());
    return out;
}

enum class SharedMode {copy, alias};
/// whether to alias or copy
extern thread_local SharedMode GlobalSharedMode;
/// vector of (shared_ptr address, respective json representation)
extern thread_local std::vector<std::pair<void const *, nlohmann::json>> GlobalSharedState;

/// RAII handle for saving/loading shared_ptrs since global objects are involved
struct Handle {
    Handle(SharedMode m) {GlobalSharedMode = m; GlobalSharedState.clear();}
    ~Handle() {GlobalSharedState.clear(); GlobalSharedMode = SharedMode::copy;}
};
/**
 * @brief Save JSON value from object that may contain shared_ptrs
 * Returns just the JSON value if json::copy or else an array [JSON value, shared_ptr JSON values]
 * @param t the object
 * @param mode either json::copy or json::alias
 */
template <class T>
nlohmann::json save_shared(T &&t, SharedMode mode=SharedMode::alias) {
    Handle handle{mode};
    nlohmann::json out = fw<T>(t);
    if (mode == SharedMode::copy) return out;
    auto shared = nlohmann::json::array();
    for (auto const &p : GlobalSharedState) shared.push_back(std::move(p.second));
    return nlohmann::json::array({std::move(out), std::move(shared)});
}
/**
 * @brief Load JSON value into an object that may contain shared_ptrs
 *
 * @param j the JSON for the object
 * @param shared the JSON for the global shared_ptr store (if null, copy each separately)
 */
template <class T>
T load_shared(nlohmann::json const &j, nlohmann::json const &shared=nullptr) {
    Handle handle{shared.is_null() ? SharedMode::copy : SharedMode::alias};
    if (!shared.is_null())
        for (auto const &j : shared) GlobalSharedState.emplace_back(nullptr, j);
    return j.get<T>();
}

}

/******************************************************************************************/

namespace nlohmann {

/**
 * @brief Overload for JSON serialization of shared_ptrs
 * If objects may contain shared_ptr and aliasing is desired (to save space)
 * use load_shared and save_shared functions above at the root. This will make shared_ptr
 * serialize to an integer, and save_shared() will return an array of length 2, with the second
 * object being the global store of shared_ptr objects indexed by that integer.
 */
template <class T>
struct adl_serializer<std::shared_ptr<T>, void> {
    static void to_json(json &j, std::shared_ptr<T> const &t) {
        auto const n = t.use_count();
        if (n == 0) {
            j = nullptr;
        } else if (nupack::GlobalSharedMode == nupack::SharedMode::copy) {
            j = *t;
        } else {
            auto &s = nupack::GlobalSharedState;
            std::size_t i;
            auto ptr = static_cast<void const *>(t.get());
            auto it = std::find_if(s.begin(), s.end(), [=](auto const &p) {return p.first == ptr;});
            if (it == s.end()) {
                i = s.size();
                s.emplace_back(ptr, *t);
            } else i = it - s.begin();
            j = static_cast<std::uint32_t>(i);
        }
    }
    static void from_json(json const &j, std::shared_ptr<T> &t) {
        if (j.is_null()) {
            t = nullptr;
        } else if (nupack::GlobalSharedMode == nupack::SharedMode::copy) {
            t = std::make_shared<T>(j.get<T>());
        } else {
            auto idx = j.get<std::uint32_t>();
            auto &l = nupack::GlobalSharedState[idx];
            if (l.first == nullptr) {
                t = std::make_shared<T>(l.second.get<T>());
                l.first = static_cast<void const *>(t.get());
            } else {
                t = *reinterpret_cast<std::shared_ptr<T> const *>(l.first);
            }
        }
    }
};

/// JSON serializer for empty object returns "{}"
template <class T>
struct adl_serializer<T, nupack::void_if<nupack::is_empty<T>>> {
    static void to_json(json &j, T const &) {j = json::object();}
    static void from_json(json const &, T &) {}
};

/// JSON Serializer for unique_ptr, gives null if it is empty
template <class T>
struct adl_serializer<std::unique_ptr<T>, void> {
    static void to_json(json &j, std::unique_ptr<T> const &value) {
        if (value) j = json::array({*value});
        else j = nullptr;
    }
    static void from_json(json const &j, std::unique_ptr<T> &value) {
        if (!j.is_null()) value = std::make_unique<T>(j.at(0).get<T>());
    }
};

/// JSON Serializer for optional, gives null if it is empty
template <class T>
struct adl_serializer<nupack::Optional<T>, void> {
    static void to_json(json &j, nupack::Optional<T> const &value) {
        if (value) j = json::array({*value});
        else j = nullptr;
    }
    static void from_json(json const &j, nupack::Optional<T> &value) {
        if (!j.is_null()) value = j.at(0).get<T>();
    }
};

/// JSON Serializer for optional, gives null if it is empty
template <class T>
struct adl_serializer<nupack::View<T>, void> {
    static void to_json(json &j, nupack::View<T> const &value) {
        j = json::array(value.begin(), value.end());
    }
    static void from_json(json const &j, nupack::View<T> &value) {
        if (j.size() != value.size()) throw std::runtime_error("incorrect array size");
        std::copy(j.begin(), j.end(), value.begin());
    }
};

/// JSON serializer for variants saves as {"index": N, "value": V}
template <class ...Ts>
struct adl_serializer<std::variant<Ts...>, void> {
    static void to_json(json &j, std::variant<Ts...> const &t) {
        j["index"] = std::uint32_t(t.index());
        nupack::fork(t, [&j](auto const &v) {j["value"] = v;});
    }
    static void from_json(json const &j, std::variant<Ts...> &t) {
        int n = j["index"].get<std::uint32_t>();
        nupack::pack<Ts...>::while_each([&](auto x) {
            if (n--) return true;
            t = j["value"].get<decltype(*x)>();
            return false;
        });
        if (n != -1) throw std::out_of_range("JSON variant index out of bounds");
    }
};

/// JSON serializer for maps where the key type is not convertible from string
template <class T>
struct adl_serializer<T, nupack::void_if<(nupack::is_map<T> && !nupack::can_convert<typename T::key_type, std::string>)>> {
    static void to_json(json &j, T const &value) {
        for (auto const &p : value) j.push_back({p.first, p.second});
    }
    static void from_json(json const &j, T &value) {
        for (auto const &p : j) value.emplace_hint(value.end(), p.at(0).template get<typename T::key_type>(),
                                                                p.at(1).template get<typename T::mapped_type>());
    }
};

/**
 * @brief Overloads for JSON treating a class as a JSON key-value dict
 *
 * @tparam T satisfies
 *      names_of(T const &) returns compile time list of names
 *      members_of(T &) returns compile time references to all members of the class
 */
template <class T>
struct adl_serializer<T, nupack::void_if<(nupack::is_nupack<T> && nupack::has_members<T> && !nupack::has_save_repr<T>)>> {
    static void to_json(json &j, T const &t) {
        j = json::object();
        nupack::for_each_zip(nupack::names_of(t), nupack::members_of(t), [&](auto n, auto const &m) {j[n] = m;});
    }
    static void from_json(json const &j, T &t) {
        nupack::for_each_zip(nupack::names_of(t), nupack::members_of(t), [&](auto n, auto &m) {j.at(n).get_to(m);});
    }
};

/**
 * @brief Overloads for JSON treating a class as a JSON value
 *
 * @tparam T satisfies
 *      declval<T const &>().save_repr() -> O such that O is convertible to JSON value
 *      declval<T &>().load_repr(O) loads the object
 */
template <class T>
struct adl_serializer<T, nupack::void_if<(nupack::is_nupack<T> && nupack::has_save_repr<T> && !nupack::has_repr_names<T>)>> {
    static void to_json(json &j, T const &t) {j = t.save_repr();}
    static void from_json(json const &j, T &t) {t.load_repr(j.get<decltype(t.save_repr())>());}
};

/**
 * @brief Overloads for JSON treating a class as a JSON key-value dict
 *
 * @tparam T satisfies
 *      declval<T const &>().repr_names() returns std::get-able list of names
 *      declval<T const &>().save_repr() returns tuple of objects (each non-reference or const-reference)
 *      declval<T &>().load_repr(...) loads those same objects objects (each non-reference or const-reference)
 */
template <class T>
struct adl_serializer<T, nupack::void_if<(nupack::is_nupack<T> && nupack::has_repr_names<T>)>> {
    static void to_json(json &j, T const &t) {
        nupack::for_each_zip(t.repr_names(), t.save_repr(), [&](auto n, auto const &m) {j[n] = m;});
    }
    template <class N, std::size_t ...Is>
    static void from_json_helper(json const &j, T &t, N const &names, nupack::indices_t<Is...>) {
        using O = decltype(t.save_repr());
        t.load_repr((j[std::get<Is>(names)].template get<nupack::no_qual<std::tuple_element_t<Is, O>>>())...);
    }
    static void from_json(json const &j, T &t) {
        from_json_helper(j, t, T::repr_names(), nupack::indices_in(T::repr_names()));
    }
};

/******************************************************************************************/

}

