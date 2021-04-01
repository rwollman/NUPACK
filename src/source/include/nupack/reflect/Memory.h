/** \file Memory.h
 * @brief Basic templates for getting memory usage of classes
 */
#pragma once
#include "../iteration/Patterns.h"
#include "Reflection.h"
#include <typeindex>

namespace nupack {

namespace memory {

/******************************************************************************************/

/// is_simple marks if a class should just be measured by sizeof() or if it can be introspected
template <class T, class=void> struct is_simple : bool_t<(is_ref_wrapper<T> || is_empty<T>
    || is_enum<T> || is_scalar<T> || is_character<T> || is_ref<T> || is_dumb_ptr<T> || std::is_pod<T>::value)> {};

template <class T> struct is_simple<T, void_if<T::simple_type::value>> : True {};
template <> struct is_simple<std::type_index, void> : True {};

/******************************************************************************************/

/** Memory implementation is always asserted to exist, so this must be defined.
 * Use is_simple for the easiest solution. Include "using simple_type = True" in the class,
 * partially specialize the is_simple struct above, use members API, or define custom impl
 */

template <class T, class=void>
struct impl {
    static_assert(!is_nupack<T> || is_simple<T>::value, "You need to specialize this if not simple");
    constexpr auto operator()(T const &) const {return sizeof(if_t<is_complete<T>, T, std::true_type>);}
    void erase(T const &) const {}
};

/******************************************************************************************/

/// For const objects, prevent memory release
template <class T> struct impl<T const, void> {
    constexpr std::size_t operator()(T const &t) const {return impl<T>()(t);}
    void erase(T const &) const {}
};

/******************************************************************************************/


template <class T> struct impl<ref_member<T>> {
    constexpr std::size_t operator()(ref_member<T> const &t) const {return impl<decay<T>>()(t.get());}
    void erase(ref_member<T> &t) const {impl<decay<T>>().release(t.get());}
};

/******************************************************************************************/

/// std::string memory implementation
template <class T> struct impl<std::basic_string<T>, void> {
    constexpr std::size_t operator()(std::basic_string<T> const &t) const {return t.capacity() * sizeof(T);}
    void erase(std::basic_string<T> &t) const {std::basic_string<T> t_; t.swap(t_);}
};

/******************************************************************************************/

/// std::pair memory implementation
template <class T> struct impl<T, void_if<is_pair<T>>> {
    template <class T_> constexpr std::size_t operator()(T_ const &t) const {
        return impl<no_cv<first_type_of<T>>>()(t.first)
             + impl<no_cv<second_type_of<T>>>()(t.second);
    }

    void erase(T &t) const {
        impl<first_type_of<T>>().release(remove_const(t.first));
        impl<second_type_of<T>>().release(remove_const(t.second));
    }
};

/******************************************************************************************/

/// std::tuple memory implementation
template <class T> struct impl<T, void_if<is_tuple<T>>> {
    template <class T_, std::size_t ...Is> constexpr auto count(T_ const &t, indices_t<Is...>) const {
        return fold(plus, 0u, impl<no_cv<tuple_type<T, Is>>>()(std::get<Is>(t))...);
    }

    template <class T_> constexpr auto operator()(T_ const &t) const {return count(t, indices_in(t));}

    template <std::size_t ...Is> void erase(T &t, indices_t<Is...>) const {
        NUPACK_UNPACK(impl<tuple_type<T, Is>>().release(remove_const(std::get<Is>(t))));
    }

    void erase(T &t) const {release(t, indices_in(t));}
};

/******************************************************************************************/

/// Class with declared members API
template <class T> struct impl<T, void_if<has_members<T>>> {
    template <class M, std::size_t ...Is>
    constexpr auto count(M const &m, indices_t<Is...>) const {
        return impl<std::tuple<no_ref<tuple_type<M, Is>>...>>()(m);
    }

    template <class T_> constexpr auto operator()(T_ const &t) const {
        return count(members_of(t), indices_in(members_of(t)));
    }

    struct releaser { // only release non-const lvalues
        template <class U> void operator()(U &u) const {impl<U>().release(u);}
        template <class U> void operator()(U const &) const {}
    };

    /// If the class contained a reference member, it should have been wrapped with std::ref or std::cref
    void erase(T &t) const {for_each(members_of(t), releaser());}
};

/******************************************************************************************/

/// Get the size in memory of an object(s) as an integer
struct Measure {
    template <class ...Ts> constexpr std::size_t operator()(Ts const &...ts) const {
        std::size_t out = fold(plus, 0u, memory::impl<Ts>()(ts)...);
        return std::max(out, (sizeof(Ts) +  ...));
    }
};

/// Release the heap memory in an object
struct Release {template <class T> void operator()(T &t) const {memory::impl<T>().erase(t);}};

/******************************************************************************************/

/// sizeof functor
template <class T> constexpr std::size_t sizeof_f(T const &) {return sizeof(T);}

/******************************************************************************************/

static constexpr auto measure = memory::Measure();
static constexpr auto release = memory::Release();

}

/******************************************************************************************/

}
