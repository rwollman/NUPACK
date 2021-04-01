/**
 * @brief C++ type-erased Variable object
 * @file Variable.h
 */

#pragma once
#include "Storage.h"
#include "Signature.h"

#include <iostream>
#include <vector>
#include <type_traits>
#include <string_view>
#include <typeindex>
#include <optional>

namespace rebind {

using TargetQualifier = Qualifier;

/******************************************************************************/

template <class SourceType, TargetQualifier Q=Value, class SFINAE=void>
struct Response; // converts type T with qualifier Q into any qualified type

template <class TargetType, class SFINAE=void>
struct Request; // makes type T, a qualified type, from a Variable

template <class T>
struct Action;

/******************************************************************************/

class Variable : protected VariableData {
    Variable(void *p, TypeIndex idx, ActionFunction act, bool s) noexcept
        : VariableData(p ? idx : TypeIndex(), p ? act : nullptr, p && s)
        {if (p) reinterpret_cast<void *&>(buff) = p;}

    Variable(Variable const &v, bool move) : VariableData(v) {
        if (v.has_value())
            act(move ? ActionType::move : ActionType::copy, v.pointer(), this);
        idx.set_qualifier(Value);
    }

    Variable request_var(Dispatch &msg, TypeIndex const &, Qualifier source) const;

public:

    Qualifier qualifier() const {return idx.qualifier();}
    void const * data() const {return pointer();}

    TypeIndex type() const {return idx;}
    ActionFunction action() const {return act;}
    bool is_stack_type() const {return stack;}

    /**************************************************************************/

    constexpr Variable() noexcept = default;

    /// Reference type
    template <class T, std::enable_if_t<!(std::is_same_v<std::decay_t<T>, T>), int> = 0>
    Variable(Type<T> t, typename SameType<T>::type reference) noexcept
        : VariableData(t, Action<std::decay_t<T>>::apply, UseStack<unqualified<T>>::value) {
            reinterpret_cast<std::remove_reference_t<T> *&>(buff) = std::addressof(reference);
        }

    /// Non-Reference type
    template <class T, class ...Ts, std::enable_if_t<(std::is_same_v<std::decay_t<T>, T>), int> = 0>
    Variable(Type<T> t, Ts &&...ts) : VariableData(t, Action<T>::apply, UseStack<T>::value) {
        static_assert(!std::is_same_v<unqualified<T>, Variable>);
        if constexpr(UseStack<T>::value) ::new (&buff) T{static_cast<Ts &&>(ts)...};
        else reinterpret_cast<T *&>(buff) = ::new T{static_cast<Ts &&>(ts)...};
    }

    template <class T, std::enable_if_t<!(std::is_same_v<std::decay_t<T>, T>), int> = 0>
    std::remove_reference_t<T> *emplace(Type<T> t, typename SameType<T>::type reference) {
        using U = unqualified<T>;
        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
        static_cast<VariableData &>(*this) = {t, Action<U>::apply, UseStack<U>::value};
        return reinterpret_cast<std::remove_reference_t<T> *&>(buff) = std::addressof(reference);
    }

    template <class T, class ...Ts, std::enable_if_t<(std::is_same_v<T, std::decay_t<T>>), int> = 0>
    T * emplace(Type<T> t, Ts &&...ts) {
        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
        static_cast<VariableData &>(*this) = {t, Action<T>::apply, UseStack<T>::value};
        if constexpr(UseStack<T>::value) return ::new (&buff) T{static_cast<Ts &&>(ts)...};
        else return reinterpret_cast<T *&>(buff) = ::new T{static_cast<Ts &&>(ts)...};
    }

    template <class T, std::enable_if_t<!std::is_base_of_v<VariableData, unqualified<T>>, int> = 0>
    Variable(T &&t) : Variable(Type<std::decay_t<T>>(), static_cast<T &&>(t)) {
        static_assert(!std::is_same_v<unqualified<T>, Variable>);
    }

    /// Take variables and reset the old ones
    // If RHS is Reference, RHS is left unchanged
    // If RHS is Value and held in stack, RHS is moved from
    // If RHS is Value and not held in stack, RHS is reset
    Variable(Variable &&v) noexcept : VariableData(static_cast<VariableData const &>(v)) {
        if (auto p = v.handle()) {
            if (stack) act(ActionType::move, p, this);
            else v.reset_data();
        }
    }

    /// Only call variable copy constructor if its lifetime is being managed
    Variable(Variable const &v) : VariableData(static_cast<VariableData const &>(v)) {
        if (auto p = v.handle()) act(ActionType::copy, p, this);
    }

    template <class T, std::enable_if_t<!std::is_base_of_v<VariableData, unqualified<T>>, int> = 0>
    Variable & operator=(T &&t) {emplace(Type<std::decay_t<T>>(), static_cast<T &&>(t)); return *this;}

    /// Only call variable move constructor if its lifetime is being managed inside the buffer
    Variable & operator=(Variable &&v) noexcept {
        // DUMP("move assign ", type(), v.type());
        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
        static_cast<VariableData &>(*this) = v;
        if (auto p = v.handle()) {
            if (stack) act(ActionType::move, p, this);
            else v.reset_data();
        }
        return *this;
    }

    Variable & operator=(Variable const &v) {
        // DUMP("copy assign ", type(), v.type());
        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
        static_cast<VariableData &>(*this) = v;
        if (auto p = v.handle()) v.act(ActionType::copy, p, this);
        return *this;
    }

    ~Variable() {
        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
    }

    /**************************************************************************/

    void reset() {
        // DUMP("reset", type());

        if (auto p = handle()) act(ActionType::destroy, p, nullptr);
        reset_data();
    }

    void assign(Variable v);

    constexpr bool has_value() const {return act;}
    explicit constexpr operator bool() const {return act;}

    /**************************************************************************/

    Variable copy() && {return {*this, qualifier() == Value || qualifier() == Rvalue};}
    Variable copy() const & {return {*this, qualifier() == Rvalue};}

    Variable reference() & {return {pointer(), idx.add(Lvalue), act, stack};}
    Variable reference() const & {return {pointer(), idx.add(Const), act, stack};}
    Variable reference() && {return {pointer(), idx.add(Rvalue), act, stack};}

    Variable request_variable(Dispatch &msg, TypeIndex const &t) const & {return request_var(msg, t, add(qualifier(), Const));}
    Variable request_variable(Dispatch &msg, TypeIndex const &t) & {return request_var(msg, t, add(qualifier(), Lvalue));}
    Variable request_variable(Dispatch &msg, TypeIndex const &t) && {return request_var(msg, t, add(qualifier(), Rvalue));}

    bool move_if_lvalue() {return idx.qualifier() == Lvalue ? idx.set_qualifier(Rvalue), true : false;}

    /**************************************************************************/

    // request reference T by custom conversions
    template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
    std::remove_reference_t<T> *request(Dispatch &msg, Type<T> t={}) const {
        DUMP("Variable.request() ", typeid(Type<T>).name(), qualifier(), " from variable ", idx);
        if (idx.matches<T>()) return target<T>();
        auto v = request_variable(msg, type_index<T>());
        if (auto p = v.template target<T>()) {msg.source.clear(); return p;}
        if (auto p = Request<T>()(*this, msg)) {msg.source.clear(); return p;}
        return nullptr;
    }

    template <class T, std::enable_if_t<!std::is_reference_v<T>, int> = 0>
    std::optional<T> request(Dispatch &msg, Type<T> t={}) const {
        static_assert(std::is_same_v<T, unqualified<T>>);
        if constexpr(std::is_same_v<T, Variable>) return *this;

        std::optional<T> out;
        if (auto p = target<T const &>()) out.emplace(*p);
        else {
            auto v = request_variable(msg, typeid(T));
            if (auto p = std::move(v).target<T &&>()) {msg.source.clear(); out.emplace(std::move(*p));}
            else if ((out = Request<T>()(*this, msg))) msg.source.clear();
        }
        // DUMP(type(), p, &buff, reinterpret_cast<void * const &>(buff), stack, typeid(p).name(), typeid(Type<T>).name());

        return out;
    }

    /**************************************************************************/

    void cast(Dispatch &msg, Type<void> t={}) const {}
    void cast(Type<void> t={}) const {}

    template <class T>
    T cast(Dispatch &msg, Type<T> t={}) const {
        if (auto p = request(msg, t)) return static_cast<T>(*p);
        throw std::move(msg).exception();
    }

    // request non-reference T by custom conversions
    template <class T>
    std::optional<T> request(Type<T> t={}) const {Dispatch msg; return request(msg, t);}

    // request non-reference T by custom conversions
    template <class T>
    T cast(Type<T> t={}) const {
        Dispatch msg;
        if (auto p = request(msg, t))
            return msg.storage.empty() ? static_cast<T>(*p) : throw std::runtime_error("contains temporaries");
        return cast(msg, t);
    }

    template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
    std::remove_reference_t<T> *target(Type<T> t={}) && {
        // DUMP(name(), typeid(Type<T>).name(), qual, stack);
        return target_pointer(t, add(qualifier(), Rvalue));
    }

    // return pointer to target if it is trivially convertible to requested type
    template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
    std::remove_reference_t<T> *target(Type<T> t={}) const & {
        // DUMP(name(), typeid(Type<T>).name(), qual, stack);
        return target_pointer(t, add(qualifier(), Const));
    }

    // return pointer to target if it is trivially convertible to requested type
    template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
    std::remove_reference_t<T> *target(Type<T> t={}) & {
        // DUMP(name(), typeid(Type<T>).name(), qual, stack);
        return target_pointer(t, add(qualifier(), Lvalue));
    }
};

/******************************************************************************/

void set_source(Dispatch &, std::type_info const &, Variable &&v);

template <TargetQualifier Q, class T>
bool get_response(Variable &out, TypeIndex const &target_type, T &&t) {
    DUMP("get_response: trying to get ", target_type, " from ", type_index<T &&>());
    if (target_type.matches<T>()) {
        DUMP("get_response: requested type matches held type");
        if constexpr(std::is_convertible_v<T &&, qualified<unqualified<T>, Q>>)
            return out = {Type<qualified<unqualified<T>, Q>>(), static_cast<T &&>(t)}, true;
    } else {
        using R = Response<unqualified<T>, Q>;
        DUMP("get_response: calling Response specialization ", type_index<R>());
        // if constexpr(std::is_invocable_v<R, Variable &, TypeIndex &&, T &&>) {
            // DUMP("get_response Response works", type_index<R>());
        bool ok = R()(out, target_type, static_cast<T &&>(t));
        DUMP("get_response: got result of type ", out.type());
        return ok;
        // }
    }
    return false;
}

/// Set out to a new variable with given qualifier dest and type idx from type T
template <class T>
bool get_response(Variable &out, TypeIndex const &target_type, T &&t) {
    DUMP("get_response: ", target_type, type_index<T &&>());
    // Switch on the runtime value of the qualifier to hit compile-time overloads
    switch (target_type.qualifier()) {
        case (Value): return get_response<Value>(out, target_type, static_cast<T &&>(t));
        case (Const): return get_response<Const>(out, target_type, static_cast<T &&>(t));
        case (Lvalue): return get_response<Lvalue>(out, target_type, static_cast<T &&>(t));
        case (Rvalue): return get_response<Rvalue>(out, target_type, static_cast<T &&>(t));
    }
    return false;
}

template <class T>
bool get_response(TypeIndex const &target_type, T &&t) {
    Variable out;
    return get_response(out, target_type, static_cast<T &&>(t));
}

/******************************************************************************/


template <class T>
struct Action {
    static_assert(std::is_same_v<unqualified<T>, T>);

    static void response(Variable &v, void *p, RequestData &&r) {
        bool ok = false;
        Dispatch &msg = *r.msg; // r is aliasing v, so save a copy of the reference
        if (r.source == Const)
            ok = get_response(v, std::move(r.type), *static_cast<T const *>(p));
        else if (r.source == Lvalue)
            ok = get_response(v, std::move(r.type), *static_cast<T *>(p));
        else if (r.source == Rvalue)
            ok = get_response(v, std::move(r.type), static_cast<T &&>(*static_cast<T *>(p)));
        else throw std::invalid_argument("source qualifier should not be Value");
        if (!ok) {
            set_source(msg, typeid(T), std::move(v)); v.reset();
        }
    }

    static void apply(ActionType a, void *p, VariableData *v) {
        if (a == ActionType::destroy) { // Delete the object (known to be non-reference)
            if constexpr(UseStack<T>::value) static_cast<T *>(p)->~T();
            else delete static_cast<T *>(p);

        } else if (a == ActionType::copy) { // Copy-Construct the object
            DUMP(v->stack, UseStack<T>::value);
            if constexpr(UseStack<T>::value) ::new(static_cast<void *>(&v->buff)) T{*static_cast<T const *>(p)};
            else reinterpret_cast<void *&>(v->buff) = ::new T{*static_cast<T const *>(p)};

        } else if (a == ActionType::move) { // Move-Construct the object (known to be on stack)
            DUMP(v->stack, UseStack<T>::value);
            if constexpr(UseStack<T>::value) // this is always known, but eliminates compile warnings
                ::new(static_cast<void *>(&v->buff)) T{std::move(*static_cast<T *>(p))};

        } else if (a == ActionType::response) { // Respond to a given type_index
            response(reinterpret_cast<Variable &>(*v), p, std::move(reinterpret_cast<RequestData &>(v->buff)));

        } else if (a == ActionType::assign) { // Assign from another variable
            // DUMP("assign", v->idx.name(), typeid(T).name(), v->qual);
            if constexpr(std::is_move_assignable_v<T>) {
                if (auto r = reinterpret_cast<Variable &&>(*v).request<T>()) {
                    // DUMP("got the assignable", v->idx.name(), typeid(T).name(), v->qual, typeid(T).name());
                    *static_cast<T *>(p) = std::move(*r);
                    reinterpret_cast<Variable &>(*v).reset(); // signifies that assignment took place
                }
            }
        }
    }
};

/******************************************************************************/

}
