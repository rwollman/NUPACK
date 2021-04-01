
#pragma once
#define BOOST_NO_IOSTREAM

#include "Memory.h"

#include "../iteration/Patterns.h"
#include "../iteration/View.h"
#include "../iteration/Search.h"
#include "../algorithms/Utility.h"
#include "../common/Runtime.h"
#include "../common/Threading.h"
#include "../common/Config.h"
#include "../common/Config.h"

#include <iostream>
#include <iomanip>
#include <mutex>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

/******************************************************************************************/

#define NUPACK_PRINT_IF(expr) template <class T> struct Printer<T, void_if<(expr)>>

namespace nupack {

/******************************************************************************************/

NUPACK_DEFINE_TYPE(is_ios_flag, decay<decltype(std::setprecision(10))>);
NUPACK_EXTEND_TYPE(is_ios_flag, decay<decltype(std::setbase(10))>);
NUPACK_EXTEND_TYPE(is_ios_flag, decay<decltype(std::setfill(' '))>);
NUPACK_EXTEND_TYPE(is_ios_flag, decay<decltype(std::setw(10))>);
NUPACK_DEFINE_VARIADIC(is_string, std::basic_string, class);

/******************************************************************************************/

namespace io {

/// Default out stream used in NUPACK (initialized to std::cout)
extern std::ostream *default_out;
/// Mutex protecting the default out stream
extern std::mutex default_out_guard;
/// Return a reference to the default out stream
inline std::ostream & out() {return *default_out;}
/// Set the default out stream (lock is used)
inline void set_out(std::ostream &os) {with_lock(default_out_guard, [&]{default_out = &os;});}

/******************************************************************************************/

struct Indent {
    int size = 4, shift = 0;
    Indent operator+() const {return {size, shift + 1};}
    friend std::ostream & operator<<(std::ostream &os, Indent t) {return os << string(t.shift * t.size, ' ');}
};

/******************************************************************************************/

template <class T, class=void> struct PrintAsContainer : False {};

struct PrintAsList : True {
    static constexpr auto begin() {return '[';}
    static constexpr auto end() {return ']';}
};

struct PrintAsSet : True {
    static constexpr auto begin() {return '{';}
    static constexpr auto end() {return '}';}
};

template <class T> struct PrintAsContainer<T, void_if<is_string<T>>> : PrintAsList {};

/******************************************************************************************/

template <class T, class=void> struct SingleLine : True {
    constexpr bool operator()(T const &) const {return true;}
};

NUPACK_UNARY_FUNCTOR(is_single_line, bool(SingleLine<std::decay_t<T>>()(t)));

/******************************************************************************************/

template <class T>
struct SingleLine<T, void_if<has_members<T> && !PrintAsContainer<T>::value && !is_streamable<T>>> {
    template <class M, std::size_t ...Is>
    static constexpr bool check(M const &m, indices_t<Is...>) {return (is_single_line(std::get<Is>(m)) && ...);}

    constexpr bool operator()(T const &t) const {
        return check(members_of(t), member_indices(t)) && memory::measure(t) < 512;
    }
};

template <class T>
struct SingleLine<T, void_if<has_len<T> && PrintAsContainer<T>::value && !is_streamable<T> && !is_character<value_type_of<T>>>> {
    bool operator()(T const &t) const {
        return !len(t) || (all_of(t, is_single_line) && memory::measure(t) < 1024);
    }
};

template <class T>
struct SingleLine<T, void_if<is_tuple<T>>> {
    template <std::size_t ...Is>
    static constexpr bool impl(T const &t, indices_t<Is...>) {
        return !sizeof...(Is) || fold(logical_or, false, is_single_line(at_c<Is>(t))...);
    }
    constexpr bool operator()(T const &t) const {return impl(t, indices_in(t));}
};

/******************************************************************************************/

template <class T, NUPACK_IF(memory::is_simple<T>::value)>
void warn_print() {}

template <class T, NUPACK_IF(!memory::is_simple<T>::value)>
[[deprecated]] void warn_print() {} // no print found

/******************************************************************************************/

template <class T, class=void>
struct Printer {
     void operator()(std::ostream &os, T const &, Indent) const {os << ~TypeName<T>() << "()";}
};

NUPACK_PRINT_IF(is_pair<T>);
NUPACK_PRINT_IF(is_tuple<T>);
NUPACK_PRINT_IF(traits::is_ref_wrapper<T>);
NUPACK_PRINT_IF(traits::is_enum<T> && !is_streamable<T>);
NUPACK_PRINT_IF(PrintAsContainer<T>::value && !is_streamable<T> && !is_character<value_type_of<T>>);
NUPACK_PRINT_IF(PrintAsContainer<T>::value && !is_streamable<T> && is_character<value_type_of<T>>);
NUPACK_PRINT_IF(has_members<T> && !PrintAsContainer<T>::value && !is_streamable<T>);
NUPACK_PRINT_IF(is_dumb_ptr<T> && !is_streamable<T> && !PrintAsContainer<T>::value);

/******************************************************************************************/

NUPACK_PRINT_IF(is_streamable<T>) {
    using uses_shift = True;
    void operator()(std::ostream &os, T const &t, Indent) const {
        if (!is_single_line(t)) os << '\n'; os << t;
    }
};

/******************************************************************************************/

NUPACK_PRINT_IF(is_pair<T>) {
    using F = decay<first_type_of<T>>;
    using S = decay<second_type_of<T>>;

    void operator()(std::ostream &os, T const &t, Indent id) const {
        os << '('; Printer<F>()(os, t.first, id); os << " : "; Printer<S>()(os, t.second, id); os << ')';
    }
};

/******************************************************************************************/

NUPACK_PRINT_IF(is_tuple<T>) {
    template <std::size_t I>
    using VT = no_ref<decltype(std::get<I>(declval<T>()))>;

    template <std::size_t ...Is>
    void helper(std::ostream &os, T const &t, Indent id, indices_t<Is...>) const {
        int n = 0;
        NUPACK_UNPACK((n++ ? os << ", " : os << '('), Printer<VT<Is>>()(os, std::get<Is>(t), +id));
        os << ')';
    }

    void operator()(std::ostream &os, T const &t, Indent id) const {helper(os, t, id, indices_in(t));}
};

/******************************************************************************************/

NUPACK_PRINT_IF(traits::is_ref_wrapper<T>) {
    void operator()(std::ostream &os, T const &t, Indent id) const {Printer<decay<decltype(t.get())>>()(os, t.get(), id);}
};

/******************************************************************************************/

NUPACK_PRINT_IF(traits::is_enum<T> && !is_streamable<T>) {
    using I = std::underlying_type_t<T>;
    void operator()(std::ostream &os, T const &t, Indent id) const {Printer<I>()(os, static_cast<I>(t), id);}
};

/******************************************************************************************/

template <> struct Printer<std::type_index> {
    void operator()(std::ostream &os, std::type_index const &t, Indent) const {os << demangle(t.name());}
};

/******************************************************************************************/

template <class T> struct Printer<ref_member<T>> {
    void operator()(std::ostream &os, ref_member<T> const &t, Indent id) const {Printer<decay<T>>()(os, t.get(), id);}
};

/******************************************************************************************/

namespace detail {
    template <class T, NUPACK_IF(has_name<T>)>
    void print_first(std::ostream &os, T const &t) {os << t.name();}

    template <class T, NUPACK_IF(!has_name<T>)>
    void print_first(std::ostream &os, T const &) {os << ~TypeName<T>();}

    template <class K, class V>
    void print_line(std::ostream &os, Indent id, K const &k, V const &v) {
        os << '\n' << id << k << ": "; Printer<V>()(os, v, id);
    }

    template <bool First, class K, class V>
    void print_part(std::ostream &os, Indent id, K const &k, V const &v) {
        (First ? os : os << ", ") << k << ": "; Printer<V>()(os, v, id);
    }
}

NUPACK_PRINT_IF(has_members<T> && !PrintAsContainer<T>::value && !is_streamable<T>) {
    template <std::size_t ...Is>
    void print_lines(std::ostream &os, T const &t, Indent id, indices_t<Is...>) const {
        NUPACK_UNPACK(detail::print_line(os, id, std::get<Is>(get_names<T>()), std::get<Is>(members_of(t))));
    }

    template <std::size_t ...Is>
    void print_parts(std::ostream &os, T const &t, Indent id, indices_t<Is...>) const {
        NUPACK_UNPACK(detail::print_part<!Is>(os, id, std::get<Is>(get_names<T>()), std::get<Is>(members_of(t))));
    }

    void operator()(std::ostream &os, T const &t, Indent id) const {
        detail::print_first(os, t);
        if (is_single_line(t)) {
            os << " {";
            print_parts(os, t, +id, member_indices(t));
            os << '}';
        } else {
            print_lines(os, t, +id, member_indices(t));
        }
    }
};

/******************************************************************************************/

NUPACK_PRINT_IF(PrintAsContainer<T>::value && !is_streamable<T> && !is_character<value_type_of<T>>) {
    void operator()(std::ostream &os, T const &t, Indent id) const {
        using VT = decay<value_type_of<T>>;
        auto const simple = is_single_line(t);
        os << PrintAsContainer<T>::begin();
        auto it = begin_of(t);
        if (it != end_of(t)) {Printer<VT>()(os, *it, +id); ++it;}
        for (; it != end_of(t); ++it) Printer<VT>()(simple ? os << ", " : os << ",\n" << +id, *it, +id);
        os << PrintAsContainer<T>::end();
    }
};

/******************************************************************************************/

NUPACK_PRINT_IF(PrintAsContainer<T>::value && !is_streamable<T> && is_character<value_type_of<T>>) {
    void operator()(std::ostream &os, T const &t, Indent) const {
        for (auto &&i : t) Printer<decay<value_type_of<T>>>()(os, i, {});
    }
};

/******************************************************************************************/

template <> struct Printer<void *> {void operator()(std::ostream &os, void *t, Indent) const {os << t;}};
template <> struct Printer<void const *> {void operator()(std::ostream &os, void const *t, Indent) const {os << t;}};

NUPACK_PRINT_IF(is_dumb_ptr<T> && !is_streamable<T> && !PrintAsContainer<T>::value) {
    void operator()(std::ostream &os, T const &t, Indent) const {os << static_cast<void const *>(&(*t));}
};

/******************************************************************************************/

template <char C=' '>
struct character {void operator()(std::ostream &os) const {os << C;}};
struct endl      {void operator()(std::ostream &os) const {os << std::endl;}};
struct new_line : character<'\n'> {};
struct tab : character<'\t'> {};
struct comma : character<','> {};
struct no_space : NoOp {};

template <class T, class=void> struct ArgPrinter {
    template <class Delim> static bool call(bool d, std::ostream &os, T const &t, Indent id) {
        if (d) Delim()(os); Printer<T>()(os, t, id); return true;
    }
};

template <> struct ArgPrinter<no_space> {
    template <class> static bool call(bool, std::ostream &os, no_space, Indent) {return false;}
};

template <class T> struct ArgPrinter<T, void_if<is_ios_flag<T>>> {
    template <class> static bool call(bool, std::ostream &os, T const &t, Indent) {os << t; return false;}
};

}

/******************************************************************************************/

/// Print to a stream with a given Indent, a given delimiter, and a given ending string
template <class Delim=io::character<' '>, class Stop=io::new_line, class ...Ts>
void print_os(io::Indent id, std::ostream &os, Ts const &...ts) {
    bool d = false;
    NUPACK_UNPACK((d = io::ArgPrinter<Ts>::template call<Delim>(d, os, ts, id)));
    Stop()(os);
}
/// Print to a stream with no Indent, a given delimiter, and a given ending string
template <class Delim=io::character<' '>, class Stop=io::new_line, class ...Ts>
void print_os(std::ostream &os, Ts const &...ts) {print_os<Delim, Stop>(io::Indent(), os, ts...);}
/// Print to a stream with no Indent, no delimiter, no ending string
template <class ...Ts>
void dump_os(std::ostream &os, Ts const &...ts) {print_os<io::no_space, io::no_space>(os, ts...);};
/// Print to a stream with an Indent, no delimiter, no ending string
template <class ...Ts>
void dump_os(io::Indent id, std::ostream &os, Ts const &...ts) {print_os<io::no_space, io::no_space>(id, os, ts...);};

/******************************************************************************************/

/// Print to global ostream, using ' ' by default as delimiter, flushing with std::endl
template <class Delim=io::character<' '>, class ...Ts>
auto print(Ts const &...ts) {
    with_lock(io::default_out_guard, [&] {
        std::ios::fmtflags f(io::out().flags());
        io::out() << std::boolalpha << std::setprecision(10);
        print_os<Delim, io::endl>(io::out(), ts...);
        io::out().flags(std::move(f));
    });
}

/// Print using new line separator
template <class ...Ts> void print_lns(Ts &&...ts) {print<io::endl>(fw<Ts>(ts)...);}

/// Expression "Echo() | arg" will print arg and return arg
struct Echo {
    template <class ...Ts> void operator()(Ts const &...ts) const {print(ts...);}
    template <class T> friend decltype(auto) operator|(T &&t, Echo const &) {print(t); return fw<T>(t);};
    template <class T> friend decltype(auto) operator|(Echo const &, T &&t) {print(t); return fw<T>(t);};
};

static constexpr auto const echo = Echo();
/// print(x, no_space, y) will print x and y with no delimiter between them
static constexpr auto const no_space = io::no_space{};

/******************************************************************************************/

inline string quoted(char const *s, bool single=false) {
    string ret;
    auto const n = std::strlen(s);
    ret.reserve(n+2);
    ret.push_back(single ? '\'' : '"');
    ret.insert(1, s, n);
    ret.push_back(single ? '\'' : '"');
    return ret;
}

inline string quoted(string s, bool single=false) {
    s.reserve(s.size() + 2);
    s.insert(0, 1, single ? '\'' : '"');
    s.push_back(single ? '\'' : '"');
    return s;
}

template <class V> string delimited_string(V const &v, string delimiter=" ") {
    std::stringstream ss; bool done = false;
    for_each(v, [&](auto const &i) {if (std::exchange(done, true)) ss << delimiter; dump_os(ss, i);});
    return ss.str();
}

template <class D, class ...Ts> string joined_string(D const &delimiter, Ts const &...ts) {
    std::stringstream ss; bool done = false;
    NUPACK_UNPACK(std::exchange(done, true) ? ss << delimiter : ss, dump_os(ss, ts));
    return ss.str();
}

/******************************************************************************************/

template <class T> struct io::PrintAsContainer<T, void_if<is_view<T>>> : PrintAsList {};

/******************************************************************************************/

namespace detail {
    // Identity operator modified for void argument
    template <class T> T const &mirror(T const &t) {return t;}
    /// Identity operator modified for void argument
    inline constexpr char const *mirror() {return "this is for the NUPACK_APPEND macro";}

    inline string stringy(char const *s) {return quoted(s);} // to print ' "x" ' rather than ' x = "x"  '
    template <class T> string stringy(T const &t) {return "do not match this";}
    inline string stringy() {return "do not match this";}
}

/******************************************************************************************/

}

/******************************************************************************************/

/// Add an argument and its value to a message. if the argument is a string just print its value
#define NUPACK_APPEND(discarded_index, discarded_data, x) \
    if (::std::string(BOOST_PP_STRINGIZE(x)) != "") { \
        if (::std::string(BOOST_PP_STRINGIZE(x)) == ::nupack::detail::stringy(x)) { \
            NUPACK_BUFFER << ::nupack::io::Indent{2, 0} << BOOST_PP_STRINGIZE(x) << '\n';} \
        else { \
            NUPACK_BUFFER << ::nupack::io::Indent{2, 0} << BOOST_PP_STRINGIZE(x) << " = "; \
            ::nupack::dump_os(NUPACK_BUFFER, ::nupack::detail::mirror(x), '\n'); \
        } \
    }

/// Capture variable or expression and print its value
#define NUPACK_CAPTURE(x) ::nupack::print("\t", ::nupack::no_space, BOOST_PP_STRINGIZE(x), " = ", x);

/// Log a message with its file location, and a variadic number of arguments
#define NUPACK_LOG(msg, ...) ({ \
    std::stringstream NUPACK_BUFFER; NUPACK_BUFFER << std::boolalpha; \
    NUPACK_BUFFER << ::nupack::message_string(NUPACK_FILE, __LINE__, msg); \
    BOOST_PP_SEQ_FOR_EACH(NUPACK_APPEND, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    ::nupack::print(NUPACK_BUFFER.str()); \
})

/// For debugging: print file location and a variadic number of arguments
#define BEEP(...) ({ \
    std::stringstream NUPACK_BUFFER; NUPACK_BUFFER << std::setprecision(16) << std::boolalpha; \
    NUPACK_BUFFER << NUPACK_FILE << ", " << __LINE__ << ": " << std::endl; \
    BOOST_PP_SEQ_FOR_EACH(NUPACK_APPEND, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)); \
    ::nupack::print(NUPACK_BUFFER.str()); \
})
