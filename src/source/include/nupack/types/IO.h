/** \file IO.h
 * @brief Contains IO functions and stream operators for many basic objects, and parameter read-in
 */
#pragma once
#include "../standard/Array.h"
#include "../standard/Vec.h"
#include "../common/Error.h"
#include "../iteration/Range.h"
#include "../iteration/Search.h"
#include "../iteration/Iterator.h"

#include <istream>
#include <string>
#include <string_view>

namespace nupack { namespace io {

bool iequals(char const *, char const *);
bool iequals(std::string const &, std::string const &);

template <class V, class S>
auto ifind(V &&v, S const &s) {
    auto it = begin_of(v);
    auto const e = end_of(v);
    for (; it != e; ++it) if (iequals(*it, s)) break;
    return it;
}

/// Returns length of deduced string minus sequence delimiters
std::size_t repeat_char_length(std::string_view s) noexcept;

/// Fills a given iterator from a string with repeated characters; returns 1 past the last iterator written to
template <class Iter>
Iter repeat_char(Iter o, std::string_view s) {
    auto b = s.begin(), e = s.end();
    if (b == e || std::isdigit(*b)) return o;
    while (b != e) {
        if (std::isdigit(*b)) {
            value_type_of<Iter> const c{*(b - 1)};
            auto const n = std::strtoul(b, const_cast<char **>(&b), 10);
            if (n > 1) o = std::fill_n(o, n - 1, c);
            else if (n == 0) --o;
        } else *(o++) = value_type_of<Iter>{*(b++)};
    }
    return o;
}

/******************************************************************************************/

/// Base class for iterator through lines of a stream
template <class S> class line_iter_base {
    S is;
    std::basic_string<typename decay<decltype(*is)>::char_type> s;
public:
    explicit line_iter_base(S p) : is(std::move(p)) {if (is) std::getline(*is, s);}
    decltype(auto) operator*() const {if (!is) throw std::out_of_range("line out of range"); return s;}
    decltype(auto) operator++() {if (is && is->good()) std::getline(*is, s); else is = nullptr; return *this;}
    bool operator==(line_iter_base const &o) const {return bool(is) == bool(o.is);}
};

/// Iterator through lines of a stream
template <class S> struct line_iter : WrapIter<line_iter<S>, line_iter_base<S>> {
    explicit line_iter(S p) : WrapIter<line_iter<S>, line_iter_base<S>>(std::move(p)) {}
};

/// Make a (Lazy) view through the lines of a stream, uses raw reference/pointer
template <class T> auto lines(std::basic_istream<T> &is) {return View<line_iter<std::basic_istream<T> *>>(&is, nullptr);}

/// Make a (Lazy) view through the lines of a stream, uses shared pointer
template <class T> auto lines(std::shared_ptr<std::basic_istream<T>> p) {
    return View<line_iter<std::shared_ptr<std::basic_istream<T>>>>(std::move(p), nullptr);
}

/// Make a (Lazy) view through the lines of a file from the name of the file
template <class T=char, class ...Ts> auto file_lines(Ts &&...ts) {
    auto p = std::make_shared<std::basic_ifstream<T>>(std::basic_ifstream<T>(fw<Ts>(ts)...));
    return lines(std::static_pointer_cast<std::istream>(std::move(p)));
}
/// Make a (Lazy) view through the lines of a string
template <class T=char, class ...Ts> auto string_lines(std::basic_string<T> const &s, Ts &&...ts) {
    auto p = std::make_shared<std::basic_stringstream<T>>(std::basic_istringstream<T>(s), fw<Ts>(ts)...);
    return lines(std::static_pointer_cast<std::istream>(std::move(p)));
}

/******************************************************************************************/

template <class V>
auto first_nonspace(V &&v) {return std::find_if_not(begin_of(v), end_of(v), is_space);}

template <class V, class F, NUPACK_IF(!(can_convert<F, value_type_of<V>>))>
bool has_content(V const &v, F &&f) {auto it = first_nonspace(v); return it != end_of(v) && f(it);}

template <class V>
bool has_content(V const &v, value_type_of<V> t) {return has_content(v, [t](auto it) {return *it != t;});}

/******************************************************************************************/

template <class T, NUPACK_IF(is_scalar<decay<T>>)>
void load_array(T &t, std::istream &is) {is >> t;}

/// Read in a std::array
template <class T, NUPACK_IF(!is_scalar<decay<T>>)>
void load_array(T &t, std::istream &is) {
    for (auto &i : t) load_array(i, is);
    if (!is.good()) NUPACK_ERROR("error in array input stream");
}

/******************************************************************************************/

template <class T> auto to_string(std::basic_istream<T> &&is) {
    std::basic_stringstream<T> buffer;
    buffer << is.rdbuf();
    return buffer.str();
}

/// Go to next line beginning with number
std::istream & go_to_number(std::istream &is);
/// Go to line after the one containing a string
std::istream & goto_line_after(std::istream &is, string const &s);
/// Skip parameter file comments
std::istream & skip_comments(std::istream &is, string const &comment_start=">");
/// Look at the next line but leave the stream in the same place
string peek(std::istream &is);
/// Check if string is on next line of stream
bool is_on_next_line(std::istream &is, string const &s);
/// Get environmental variable
string get_env(string key);

/******************************************************************************************/

/// Low level dp_to_pairs taking int * of length n, null terminated char *
void to_pairs(iseq *v, iseq const n, char const *dp);

template <class V=vec<iseq>>
V to_pairs(std::string_view dp) {
    auto n = repeat_char_length(dp);
    V v(n, {});
    if (n) to_pairs(&(v[0]), n, dp.begin());
    return v;
}

/******************************************************************************************/

/// Convert single strand complex pair array to dot-parens
template <class V> string to_dp(V const &pairs) {
    string ret(len(pairs), '_');
    for (auto i : indices(pairs)) {
        if (pairs[i] > i) ret[i] = '(';
        else if (pairs[i] < i) ret[i] = ')';
        else ret[i] = '.';
    }
    return ret;
}

/******************************************************************************************/

/// Convert multiple strand complex pair array to dot-parens
template <class V1, class V2> string to_dp(V1 const &pairs, V2 const &nicks) {
    string dp;
    auto b = begin_of(nicks), e = end_of(nicks);
    if (b != e && *b == 0) ++b; // discard beginning 0
    if (b != e && *std::prev(e) == len(pairs)) --e; // discard trailing sum
    dp.reserve(len(pairs) + std::distance(b, e));
    for (auto i : indices(pairs)) {
        if (pairs[i] > i) dp.push_back('(');
        else if (pairs[i] < i) dp.push_back(')');
        else dp.push_back('.');
        if (b != e && i + 1 == *b) {dp.push_back('+'); ++b;}
    }
    return dp;
}

/******************************************************************************************/

}}
