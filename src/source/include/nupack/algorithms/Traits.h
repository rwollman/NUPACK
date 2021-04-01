/**
 * @brief Definitions of compile-time traits
 *
 * @file Traits.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "TypeSupport.h"
#include "Macro.h"

#include <iterator>
#include <complex>
#include <string_view>

namespace nupack {

NUPACK_DEFINE_TEMPLATE(is_pair, std::pair, class, class);
NUPACK_DEFINE_TEMPLATE(is_like_tuple, std::pair, class, class);
NUPACK_DEFINE_VARIADIC(is_tuple, std::tuple, class);
NUPACK_EXTEND_VARIADIC(is_like_tuple, std::tuple, class);
NUPACK_DEFINE_VARIADIC(is_ref_wrapper, std::reference_wrapper, class);

namespace traits {

/******************************************************************************************/

/// the first type is identical to any of the latter types
template <class T, class ...Us> static constexpr bool is_same = any_of_c<std::is_same<T, Us>::value...>;
template <class T, class ...Us> static constexpr bool is_like = is_same<decay<T>, decay<Us>...>;

/// SFINAE for is_same on types
template <class T, class ...Us> using enable_same = void_if<is_same<T, Us...>>;
/// take the first type that is not void
template <class T, class U> using nonvoid = if_t<is_same<T, void>, U, T>;


/******************************************************************************************/

/// is_complex contains the underlying type as extra information
template <class T> struct is_complex : std::false_type {using type = T;};
template <class T> struct is_complex<std::complex<T>> : std::true_type {using type = T;};

template <class F, class ...Ts> std::false_type can_call_f(...);
template <class F, class ...Ts> auto can_call_f(int) -> decltype(declval<F>()(declval<Ts>()...), std::true_type());

template <class F, class T> std::false_type is_subscriptable_f(...);
template <class F, class T> auto is_subscriptable_f(int) -> decltype(declval<F>()[declval<T>()], std::true_type());

template <class F, class T> std::false_type can_assign_from_f(...);
template <class F, class T> auto can_assign_from_f(int) -> decltype(declref<F>() = declval<T>(), std::true_type());

template <class F, class ...Ts> std::false_type is_findable_f(...);
template <class F, class ...Ts> auto is_findable_f(int) -> decltype(declval<F>().find(declval<Ts>()...), std::true_type());

template <class T> std::false_type can_arrow_f(...);
template <class T> auto can_arrow_f(int) -> decltype(std::declval<T &>().operator->(), std::true_type());
template <class T> auto can_arrow_f(int) -> decltype(void_if<std::is_pointer<T>::value>(), std::true_type());

template <class T> std::false_type is_complete_f(...);
template <class T> std::true_type is_complete_f(int(*)[sizeof(T)]);

template <class T, class=void> struct is_character_t : std::false_type {};
template <class T> struct is_character_t<T, enable_same<T, char, wchar_t, char16_t, char32_t>> : std::true_type {};
template <class T> static constexpr bool is_character = is_character_t<decay<T>>::value;

template <class U> struct TypeCheck {
    template <class T> friend decltype(auto) operator|(T &&t, TypeCheck) {
        static_assert(is_same<U, T &&>, "Type check failed");
        return fw<T>(t);
    }
};

template <> struct TypeCheck<void> {
    template <class T> friend auto operator|(T &&t, TypeCheck) {
        static_assert(decltype(fw<T>(t))::checking_type, "This is the type");
    }
};

template <class ...> struct matches_signature_t : std::false_type {};
template <class F, class ...Args> struct matches_signature_t<F, std::result_of_t<F(Args...)>, Args...> : std::true_type {};
template <class ...Ts> static constexpr bool matches_signature = matches_signature_t<Ts...>::value;

template <class T, class=void> struct is_cref_wrapper_t : std::false_type {};
template <class T> struct is_cref_wrapper_t<std::reference_wrapper<T>, void_if<std::is_const<T>::value>> : std::true_type {};
template <class T> static constexpr bool is_cref_wrapper = is_cref_wrapper_t<T>::value;

/******************************************************************************************/

template <std::size_t> struct uint_of_size_t;
template <> struct uint_of_size_t<1> {using type = std::uint8_t;};
template <> struct uint_of_size_t<2> {using type = std::uint16_t;};
template <> struct uint_of_size_t<4> {using type = std::uint32_t;};
template <> struct uint_of_size_t<8> {using type = std::uint64_t;};

template <std::size_t> struct int_of_size_t;
template <> struct int_of_size_t<1> {using type = std::int8_t;};
template <> struct int_of_size_t<2> {using type = std::int16_t;};
template <> struct int_of_size_t<4> {using type = std::int32_t;};
template <> struct int_of_size_t<8> {using type = std::int64_t;};

template <std::size_t N> using int_of_size = typename int_of_size_t<N>::type;
template <std::size_t N> using uint_of_size = typename uint_of_size_t<N>::type;

/******************************************************************************************/

template <class T>
static constexpr bool is_complete = decltype(is_complete_f<T>(nullptr))::value;

template <class F, class ...Ts>
static constexpr bool can_call = decltype(can_call_f<F, Ts...>(int{}))::value;
template <class F, class T>
static constexpr bool is_subscriptable = decltype(is_subscriptable_f<F, T>(int{}))::value;
template <class F, class T>
static constexpr bool can_assign_from = decltype(can_assign_from_f<F, T>(int{}))::value;
template <class F, class ...Ts>
static constexpr bool is_findable = decltype(is_findable_f<F, Ts...>(int{}))::value;
template <class T>
static constexpr bool can_arrow = decltype(can_arrow_f<T>(int{}))::value;

/******************************************************************************************/

#define NUPACK_TMP(NAME, STD, ...) static constexpr bool NAME = std::STD<__VA_ARGS__>::value
template <class T>              NUPACK_TMP(is_scalar,         is_scalar,             decay<T>);
template <class T>              NUPACK_TMP(is_signed,         is_signed,             decay<T>);
template <class T>              NUPACK_TMP(is_unsigned,       is_unsigned,           decay<T>);
template <class T>              NUPACK_TMP(is_const,          is_const,              T);
template <class T>              NUPACK_TMP(is_enum,           is_enum,               T);
template <class T>              NUPACK_TMP(is_class,          is_class,              T);
template <class T>              NUPACK_TMP(is_arithmetic,     is_arithmetic,         decay<T>);
template <class T>              NUPACK_TMP(is_integral,       is_integral,           decay<T>);
template <class T>              NUPACK_TMP(is_floating_point, is_floating_point,     decay<T>);
template <class T>              NUPACK_TMP(is_empty,          is_empty,              decay<T>);
template <class T>              NUPACK_TMP(is_pointer,        is_pointer,            decay<T>);
template <class T>              NUPACK_TMP(is_ref,            is_reference,          T);
template <class T>              NUPACK_TMP(is_lref,           is_lvalue_reference,   T);
template <class T>              NUPACK_TMP(is_cref,           is_lvalue_reference,   T) && is_const<no_ref<T>>;
template <class T>              NUPACK_TMP(is_rref,           is_rvalue_reference,   T);
template <class T, class U>     NUPACK_TMP(is_t,              is_same,               T, U);
template <class T, class U>     NUPACK_TMP(derives_from,      is_base_of,            U, T);
template <class T, class ...Us> NUPACK_TMP(can_construct,     is_constructible,      T, Us...);
template <class T, class U>     NUPACK_TMP(can_convert,       is_convertible,        T, U);
template <class T>              NUPACK_TMP(can_copy,          is_copy_constructible, T);
template <class T>              NUPACK_TMP(can_move,          is_move_constructible, T);
template <class T>              NUPACK_TMP(can_copy_assign,   is_copy_assignable,    T);
template <class T>              NUPACK_TMP(can_move_assign,   is_move_assignable,    T);
#undef NUPACK_TMP

/******************************************************************************************/

// template <class T> static constexpr bool is_like_string = is_same<decay<T>, std::string, char *, char const *>;

/******************************************************************************************/

template <class T=void> auto type_check() {return TypeCheck<T>();}
template <class T> auto type_check(T &&t) {return TypeCheck<void>() | fw<T>(t);}

template <class T>
T what_type() {auto answer = T::this_is_the_type; return declval<T>();}

template <class T>
T what_type(T &&) {auto answer = T::this_is_the_type; return declval<T>();}

template <class ...Ts, NUPACK_IF(sizeof...(Ts) != 1)>
void what_type(Ts &&...) {auto answer = std::tuple<Ts...>::this_is_the_type;}

/******************************************************************************************/

template <class T> using remove_rref  = if_t<is_rref<T>, no_ref<T>, T>;
template <class T> using remove_lref  = if_t<is_lref<T>, no_ref<T>, T>;
template <class T> using to_reference = if_t<is_ref<T>, T, T &&>;

/******************************************************************************************/

template <class T> using iterator_category_of = typename decay<T>::iterator_category;
template <class T> using const_iterator_of =    typename decay<T>::const_iterator;
template <class T> using iterator_of =          typename decay<T>::iterator;
template <class T> using type_in =              typename decay<T>::type;
template <class T> using first_type_of =        typename decay<T>::first_type;
template <class T> using second_type_of =       typename decay<T>::second_type;
template <class T> using element_type_of =      typename decay<T>::element_type;
template <class T> using key_type_of =          typename decay<T>::key_type;
template <class T> using base_type_of =         typename decay<T>::base_type;
template <class T> using result_type_of =       typename decay<T>::result_type;
template <class T> using unsigned_type_of =     typename std::make_unsigned_t<T>;
template <class T> using signed_type_of =       typename std::make_signed_t<T>;

}

/// Should probably use std::iterator_traits for some of these...
NUPACK_TEMPLATE_FALLBACK(subtracted_type_of,   typename decay<T>::difference_type,   decltype(declval<T>() - declval<T>()));
NUPACK_TEMPLATE_FALLBACK(difference_type_of,   subtracted_type_of<T>,                std::ptrdiff_t);
NUPACK_TEMPLATE_FALLBACK(value_type_of,        typename decay<T>::value_type,        decay<decltype(*declref<T>())>);
NUPACK_TEMPLATE_FALLBACK(reference_type_of,    typename decay<T>::reference_type,    traits::to_reference<decltype(*declref<T>())>);
NUPACK_TEMPLATE_FALLBACK(size_type_of,         typename decay<T>::size_type,         std::make_unsigned_t<difference_type_of<T>>);
NUPACK_TEMPLATE_FALLBACK(pointer_type_of,      typename decay<T>::pointer,           decltype(&(*declref<T>())));

/******************************************************************************************/

NUPACK_DETECT_BOOL_OP(has_lt, can_lt, <);
NUPACK_DETECT_BOOL_OP(has_gt, can_gt, >);
NUPACK_DETECT_BOOL_OP(has_le, can_le, <=);
NUPACK_DETECT_BOOL_OP(has_ge, can_ge, >=);
NUPACK_DETECT_BOOL_OP(has_eq, can_eq, ==);
NUPACK_DETECT_BOOL_OP(has_ne, can_ne, !=);

NUPACK_DETECT(can_add,                   decltype(declval<T>() + declval<T>()));
NUPACK_DETECT(can_subtract,              decltype(declval<T>() - declval<T>()));
NUPACK_DETECT(can_multiply,              decltype(declval<T>() * declval<T>()));
NUPACK_DETECT(can_divide,                decltype(declval<T>() / declval<T>()));

NUPACK_DETECT(can_increment,             decltype(++declref<T>()));
NUPACK_DETECT(can_decrement,             decltype(--declref<T>()));
NUPACK_DETECT(can_post_increment,        decltype(declref<T>()++));
NUPACK_DETECT(can_post_decrement,        decltype(declref<T>()--));
NUPACK_DETECT(can_plus_one,              decltype(declval<T>() + 1u));
NUPACK_DETECT(can_minus_one,             decltype(declval<T>() - 1u));
NUPACK_DETECT(can_dereference,           decltype(*declref<T>()));

NUPACK_DETECT(has_value_type,            typename T::value_type);
NUPACK_DETECT(has_size_type,             typename T::size_type);
NUPACK_DETECT(has_difference_type,       typename T::difference_type);
NUPACK_DETECT(has_pointer_type,          typename T::pointer);
NUPACK_DETECT(has_reference_type,        typename T::reference);

NUPACK_DETECT(has_capacity,              decltype(declref<T>().capacity()));
NUPACK_DETECT(has_reserve,               decltype(declref<T>().reserve(1u)));
NUPACK_DETECT(has_resize,                decltype(declref<T>().resize(1u)));
NUPACK_DETECT(has_what,                  decltype(declval<T>().what()));
NUPACK_DETECT(has_begin,                 decltype(declval<T>().begin()));
NUPACK_DETECT(has_data,                  decltype(declval<T>().data()));
NUPACK_DETECT(has_end,                   decltype(declval<T>().end()));
NUPACK_DETECT(is_ptr,                    decltype(static_cast<void const *>(&(*declval<T>()))));
NUPACK_DETECT(is_fancy_ptr,              decltype(declval<T>().reset()));
NUPACK_DETECT(is_dumb_ptr,               void_if<traits::is_ptr<T> && !traits::is_fancy_ptr<T>>);
NUPACK_DETECT(is_random_access_iterator, decltype(declval<T>()[2]));
NUPACK_DETECT(has_member_swap,           decltype(declref<T>().swap(declref<T>())));
NUPACK_DETECT(has_size,                  decltype(static_cast<std::size_t>(declref<T const>().size())));

/******************************************************************************************/

template <class T, NUPACK_IF( traits::is_ptr<decay<T>>)> decltype(auto) rm_ptr(T &&t) {return *(fw<T>(t));}
template <class T, NUPACK_IF(!traits::is_ptr<decay<T>>)> decltype(auto) rm_ptr(T &&t) {return fw<T>(t);}

/******************************************************************************************/

template <class T, class U>
auto check_t(U &&u) -> void_if<traits::can_convert<decltype(u), T>>;

template <class U, class T>
decltype(auto) cast_this(T *t) {return static_cast<if_t<traits::is_const<T>, U const &, U &>>(*t);}

/******************************************************************************************/

template <class O, class T>
void variadic_assign(O &o, T &&t) {o = fw<T>(t);}

template <class O, class ...Ts, NUPACK_IF(sizeof...(Ts) != 1)>
void variadic_assign(O &o, Ts &&...ts) {o = O{fw<Ts>(ts)...};}

/******************************************************************************************/

/// this registers namespaces to nupack via ADL so we can pick up which things are ours
template <class T> std::true_type this_namespace_is_nupack(T);

#define NUPACK_NAMESPACE(NAME) namespace NAME {using ::nupack::this_namespace_is_nupack;}

NUPACK_NAMESPACE(detail);
NUPACK_NAMESPACE(traits);

using namespace traits;

/******************************************************************************************/

}

namespace nupack_adl {
    template <class T, class=void> struct is_nupack : std::false_type {};
    template <class T> struct is_nupack<T, nupack::void_if<decltype(this_namespace_is_nupack((T) nullptr))::value>> : std::true_type {};
}

namespace nupack {
    /// whether the type was declared in a nupack-registered namespace
    template <class T> static constexpr bool is_nupack = ::nupack_adl::is_nupack<decay<T> *>::value;
    static_assert(!is_nupack<int>, "namespace deduction failed");
    static_assert(is_nupack<traits::TypeCheck<void>>, "namespace deduction failed");
}
