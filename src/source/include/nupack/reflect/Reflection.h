#pragma once
#include "../algorithms/Pack.h"
#include "../iteration/Patterns.h"
#include <array>

/******************************************************************************************/

namespace nupack {

NUPACK_DETECT(has_save_repr,  decltype(declref<T const>().save_repr()));
NUPACK_DETECT(has_load_repr,  decltype(declref<T>().load_repr()));
NUPACK_DETECT(has_repr_names, decltype(T::repr_names()));

NUPACK_DETECT(has_members,    decltype(declref<T>().members()));
NUPACK_DETECT(has_names,      decltype(T::names()));
NUPACK_DETECT(has_name,       decltype(declref<T const>().name()));

NUPACK_DETECT(is_streamable,  decltype(declref<std::ostream>() << declref<no_ref<T> const>()));

/******************************************************************************************/

#define NUPACK_PUBLIC_IMPL(r, cls, i, elem) doc.method(t, "." BOOST_PP_STRINGIZE(elem), &cls::elem);
#define NUPACK_PUBLIC(cls, ...) BOOST_PP_SEQ_FOR_EACH_I(NUPACK_PUBLIC_IMPL, cls, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

/******************************************************************************************/

template <class T>
class ref_member {
    no_cv<T> *value;
public:
    using type = T;
    explicit constexpr ref_member(no_cv<T> *t) : value(t) {}
    no_cv<T> & get() const {return *value;}
    friend void swap(ref_member &self, ref_member &other) {swap(*self.value, *other.value);}
};

template <class T> auto cref_member(T &t) {return ref_member<T const>{&t};}
template <class T> auto lref_member(T &t) {return ref_member<T>{&t};}

namespace traits {
    template <class T> struct is_lref_member_t : False {};
    template <class T> struct is_lref_member_t<ref_member<T>> : bool_t<!is_const<T>> {};
    template <class T> struct is_cref_member_t : False {};
    template <class T> struct is_cref_member_t<ref_member<T>> : bool_t<is_const<T>> {};
}
template <class T> static constexpr bool is_lref_member = traits::is_lref_member_t<T>::value;
template <class T> static constexpr bool is_cref_member = traits::is_cref_member_t<T>::value;

/******************************************************************************************/

/// make_members is basically std::tie but removes rvalue references
template <class ...Cs> constexpr auto make_members(Cs &&...cs) {
    return std::tuple<remove_rref<decltype(cs)>...>{fw<Cs>(cs)...};
}

// /// extend_members appends to the members tuple
// template <class ...Ts> struct extend_members_t {
//     template <class C, class ...Cs> constexpr auto operator()(C *c, Cs &&...cs) const {
//         return std::tuple_cat(static_cast<Ts &>(*c).members()..., std::tuple<remove_rref<decltype(cs)>...>{fw<Cs>(cs)...});
//     }
// };
// template <class ...Ts> struct extend_members_t<pack<Ts...>> : extend_members_t<Ts...> {};
// template <class ...Ts> static constexpr auto extend_members = extend_members_t<Ts...>();

/******************************************************************************************/

/// make_names returns an array of char const *
template <class ...Cs> constexpr auto make_names(Cs ...cs) {
    return std::array<char const *, sizeof...(Cs)>{{cs...}};
}

// /// combine_names adds char const *s to the names array
// template <class A, class B, class ...Cs> constexpr auto combine_names(A const &a, B const &b, Cs const &...cs);
// template <class A> constexpr auto combine_names(A const &a) {return a;}

// template <class A, class B, std::size_t ...Is, std::size_t ...Js, class ...Cs>
// constexpr auto combine_names_impl(A const &a, indices_t<Is...>, B const &b, indices_t<Js...>, Cs const &...cs) {
//     return combine_names(std::array<char const *, sizeof...(Is) + sizeof...(Js)>{{at<Is>(a)..., at<Js>(b)...}}, cs...);
// }
// /// combine_names adds char const *s to the names array
// template <class A, class B, class ...Cs> constexpr auto combine_names(A const &a, B const &b, Cs const &...cs) {
//     return combine_names_impl(a, indices_in(a), b, indices_in(b), cs...);
// }

// /// extend_members appends to the members tuple
// template <class ...Ts> struct extend_names_t {
//     template <class ...Cs> constexpr auto operator()(Cs &&...cs) const {
//         return combine_names(Ts::names()..., std::array<char const *, sizeof...(Cs)>{{fw<Cs>(cs)...}});
//     }
// };
// template <class ...Ts> struct extend_names_t<pack<Ts...>> : extend_names_t<Ts...> {};
// template <class ...Ts> static constexpr auto extend_names = extend_names_t<Ts...>();

template <class T, class A, class B, std::size_t ...Is, std::size_t ...Js>
constexpr std::array<T, sizeof...(Is) + sizeof...(Js)> array_cat_impl(A &&a, B &&b, indices_t<Is...>, indices_t<Js...>) {
    return {std::get<Is>(fw<A>(a))..., std::get<Js>(fw<B>(b))...};
}

template <class T, class A, class B>
constexpr auto array_cat(A &&a, B &&b) {
    return array_cat_impl<T>(fw<A>(a), fw<B>(b), indices_in<A>(), indices_in<B>());
}

/******************************************************************************************/

NUPACK_UNARY_FUNCTOR(names_of, decay<decltype(t)>::names());

template <class T> constexpr auto member_indices() {return indices_in(decay<T>::names());}
template <class T> constexpr auto member_indices(T const &) {return member_indices<T>();}

class members_of_t {
    template <class M, std::size_t ...Is>
    static auto constify(M &&m, indices_t<Is...>) {return std::tie(add_const(at<Is>(m))...);}

    template <class M, std::size_t ...Is>
    static auto moveify(M &&m, indices_t<Is...>) {return std::forward_as_tuple(std::move(at<Is>(m))...);}

public:

    template <class T> auto operator() (T &t) const {return t.members();}

    template <class T> auto operator() (T const &t) const {
        return constify(remove_const(t).members(), member_indices<T>());
    }

    template <class T> auto operator() (T &&t) const {
        return moveify(remove_const(t).members(), member_indices<T>());
    }
};

static constexpr auto members_of = members_of_t();

/******************************************************************************************/

template <class ...Ts> static constexpr bool is_nothrow_swap = all_of_c<noexcept(swap(declref<Ts>(), declref<Ts>()))...>;

template <class T, class U, std::size_t ...Is>
void swap_all(T &&t, U &&u, indices_t<Is...>) noexcept(is_nothrow_swap<decltype(at<Is>(t))...>) {
    NUPACK_UNPACK(swap(at<Is>(t), at<Is>(u)));
}

template <class T, class U>
void swap_all(T &&t, U &&u) noexcept(noexcept(swap_all(fw<T>(t), fw<U>(u), decltype(indices_in(t))()))) {
    return swap_all(fw<T>(t), fw<U>(u), decltype(indices_in(t))());
}

template <class T>
static constexpr bool has_nothrow_swap_members = noexcept(swap_all(declval<T>().members(), declval<T>().members(), member_indices<T>()));

template <class T>
void swap_members(T &&t1, T &&t2) noexcept(has_nothrow_swap_members<T>) {
    swap_all(t1.members(), t2.members(), member_indices<T>());
}

/******************************************************************************************/

/// get_names gets the member names of a class by default
template <class U, NUPACK_IF(has_names<U>)>
static constexpr auto get_names() {return U::names();}
/// if a class has no declared names, get_names makes up its own names
template <class U, NUPACK_IF(!has_names<U>)>
static auto get_names() {
    std::array<std::string, std::tuple_size<decltype(declval<U>().members())>::value> ret;
    izip(ret, [](auto i, auto &n) {n = "unnamed_" + std::to_string(i);});
    return ret;
}

template <class T, class F>
void for_each_member(T &&t, F &&f) {
    auto mems = members_of(fw<T>(t));
    for_each_index<(std::tuple_size<decltype(mems)>::value)>([&](auto i) {return f(i, at<decltype(i)::value>(mems));});
}

/******************************************************************************************/

struct Empty {
    auto members() {return make_members();}
    static constexpr auto names() {return make_names();}
    static constexpr auto accesses() {return make_members();}

    constexpr bool operator==(Empty const &) const {return true;}
    constexpr bool operator<=(Empty const &) const {return true;}
    constexpr bool operator>=(Empty const &) const {return true;}

    constexpr bool operator!=(Empty const &) const {return false;}
    constexpr bool operator< (Empty const &) const {return false;}
    constexpr bool operator> (Empty const &) const {return false;}
};

struct MemberComparable {};
struct MemberWeaklyOrdered {};
struct MemberOrdered : MemberWeaklyOrdered, MemberComparable {};

namespace traits {
    template <class T, class=void>
    struct is_member_comparable : bool_t<derives_from<T, MemberComparable>> {};
    template <class T, class=void>
    struct is_member_weakly_ordered : bool_t<derives_from<T, MemberWeaklyOrdered>> {};
    template <class T, class=void>
    struct is_member_ordered : bool_t<derives_from<T, MemberOrdered>> {};

    template <class T> struct is_member_comparable<T, void_if<T::is_member_comparable::value>> : True {};
    template <class T> struct is_member_comparable<T, void_if<T::is_member_ordered::value>> : True {};
    template <class T> struct is_member_weakly_ordered<T, void_if<T::is_member_weakly_ordered::value>> : True {};
    template <class T> struct is_member_weakly_ordered<T, void_if<T::is_member_ordered::value>> : True {};
    template <class T> struct is_member_ordered<T, void_if<T::is_member_ordered::value>> : True {};
}

template <class T, NUPACK_IF(traits::is_member_comparable<T>::value), NUPACK_IF(has_eq<decltype(members_of(declref<T const>()))>)>
constexpr bool operator==(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? true : members_of(i) == members_of(j);}

template <class T, NUPACK_IF(traits::is_member_comparable<T>::value), NUPACK_IF(has_ne<decltype(members_of(declref<T const>()))>)>
constexpr bool operator!=(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? false : members_of(i) != members_of(j);}

template <class T, NUPACK_IF(traits::is_member_weakly_ordered<T>::value), NUPACK_IF(has_lt<decltype(members_of(declref<T const>()))>)>
constexpr bool operator<(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? false : members_of(i) < members_of(j);}

template <class T, NUPACK_IF(traits::is_member_weakly_ordered<T>::value), NUPACK_IF(has_gt<decltype(members_of(declref<T const>()))>)>
constexpr bool operator>(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? false : members_of(i) > members_of(j);}

template <class T, NUPACK_IF(traits::is_member_ordered<T>::value), NUPACK_IF(has_le<decltype(members_of(declref<T const>()))>)>
constexpr bool operator<=(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? true : members_of(i) <= members_of(j);}

template <class T, NUPACK_IF(traits::is_member_ordered<T>::value), NUPACK_IF(has_ge<decltype(members_of(declref<T const>()))>)>
constexpr bool operator>=(T const &i, T const &j) {return std::addressof(i) == std::addressof(j) ? true : members_of(i) >= members_of(j);}

/******************************************************************************************/

template <class From, class To>
struct BaseCast {
    To const & operator()(From const &f) const {return static_cast<To const &>(f);}
    To &       operator()(From &f)       const {return static_cast<To &>(f);}
    To &&      operator()(From &&f)      const {return static_cast<To &&>(f);}
};

}

#define NUPACK_MEMBERS(...) auto members() {return ::nupack::make_members(__VA_ARGS__);}

#define NUPACK_NAMES_IMPL(r, data, i, elem) BOOST_PP_COMMA_IF(i) (("." BOOST_PP_STRINGIZE(elem)) + 1)
#define NUPACK_NAMES_F(...) ::nupack::make_names(BOOST_PP_SEQ_FOR_EACH_I(NUPACK_NAMES_IMPL, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))
#define NUPACK_NAMES(...) static constexpr auto names() {return NUPACK_NAMES_F(__VA_ARGS__);}

#define NUPACK_ACCESSES_IMPL(r, cls, i, elem) BOOST_PP_COMMA_IF(i) &cls::elem
#define NUPACK_ACCESSES_F(cls, ...) ::nupack::make_members(BOOST_PP_SEQ_FOR_EACH_I(NUPACK_ACCESSES_IMPL, cls, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))
#define NUPACK_ACCESSES(cls, ...) static constexpr auto accesses() {return NUPACK_ACCESSES_F(cls, __VA_ARGS__);}

#define NUPACK_REFLECT(cls, ...) NUPACK_MEMBERS(__VA_ARGS__) NUPACK_NAMES(__VA_ARGS__) NUPACK_ACCESSES(cls, __VA_ARGS__)


/******************************************************************************************/

#define NUPACK_EXTEND_MEMBERS(base, ...) auto members() {return std::tuple_cat(base::members(), ::nupack::make_members(__VA_ARGS__));}
#define NUPACK_EXTEND_NAMES(base, ...) static constexpr auto names() {return ::nupack::array_cat<char const *>(base::names(), NUPACK_NAMES_F(__VA_ARGS__));}
#define NUPACK_EXTEND_ACCESSES(cls, base, ...) static constexpr auto accesses() {return std::tuple_cat(base::accesses(), NUPACK_ACCESSES_F(cls, __VA_ARGS__));}

#define NUPACK_EXTEND_REFLECT(cls, base, ...) NUPACK_EXTEND_MEMBERS(base, __VA_ARGS__) NUPACK_EXTEND_NAMES(base, __VA_ARGS__) NUPACK_EXTEND_ACCESSES(cls, base, __VA_ARGS__)

#define NUPACK_REFLECT_BASE(cls, base, ...) \
    NUPACK_MEMBERS(static_cast<base&>(*this), __VA_ARGS__) \
    NUPACK_NAMES(base, __VA_ARGS__) \
    static constexpr auto accesses() {return std::tuple_cat(std::make_tuple(::nupack::BaseCast<cls, base>()), NUPACK_ACCESSES_F(cls, __VA_ARGS__));}

    // static constexpr auto accesses() {return std::tuple_cat(NUPACK_ACCESSES_F(cls, ...));}
