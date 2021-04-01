/** \file IO.cc
 * @brief Defines IO functions and stream operators for many basic objects, and parameter read-in
 */

#include <nupack/types/PairList.h>
#include <nupack/types/Sequence.h>
#include <nupack/types/IO.h>
#include <nupack/types/LRU.h>
#include <nupack/types/Complex.h>

#include <nupack/reflect/Serialize.h>
#include <nupack/reflect/Print.h>
#include <nupack/iteration/Patterns.h>
#include <nupack/algorithms/Utility.h>

#include <stack>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>


namespace nupack {

constexpr std::array<char, 16> Base::names;

constexpr std::array<std::array<bool, 4>, 16> Base::masks;

constexpr Base::from_index_t const Base::from_index;

std::array<std::discrete_distribution<BaseIndex>, 15> const Base::distributions = {{
    discrete_distribution<BaseIndex>(Base::masks[0]),
    discrete_distribution<BaseIndex>(Base::masks[1]),
    discrete_distribution<BaseIndex>(Base::masks[2]),
    discrete_distribution<BaseIndex>(Base::masks[3]),
    discrete_distribution<BaseIndex>(Base::masks[4]),
    discrete_distribution<BaseIndex>(Base::masks[5]),
    discrete_distribution<BaseIndex>(Base::masks[6]),
    discrete_distribution<BaseIndex>(Base::masks[7]),
    discrete_distribution<BaseIndex>(Base::masks[8]),
    discrete_distribution<BaseIndex>(Base::masks[9]),
    discrete_distribution<BaseIndex>(Base::masks[10]),
    discrete_distribution<BaseIndex>(Base::masks[11]),
    discrete_distribution<BaseIndex>(Base::masks[12]),
    discrete_distribution<BaseIndex>(Base::masks[13]),
    discrete_distribution<BaseIndex>(Base::masks[14])
}};

/******************************************************************************************/

auto sequence_separator() {return boost::is_any_of(",+ \n\t");}

namespace io {

bool iequals(char const *s, char const *t) {return boost::iequals(s, t);}
bool iequals(std::string const &s, std::string const &t) {return boost::iequals(s, t);}

std::ostream *default_out = &std::cout;
std::mutex default_out_guard;

bool is_on_next_line(std::istream &is, string const &s) {
    string line;
    auto cur = is.tellg();
    if (!is.good()) return false;
    std::getline(is, line);
    bool ret = (line.find(s) != line.npos);
    is.seekg(cur);
    return ret;
}

/******************************************************************************************/

std::istream & go_to_number(std::istream &is) {
    auto cur = is.tellg();
    string line;
    while (std::getline(is, line)) {
        auto number_pos = line.find_first_of("+-.1234567890");
        if (number_pos != line.npos)
            if (line.find_first_not_of(" \n\r\t") == number_pos) return is.seekg(cur);
        cur = is.tellg();
    }
    throw std::runtime_error("Reached end of file while looking for numbers");
}

/******************************************************************************************/

std::istream & goto_line_after(std::istream &is, string const &s){
    NUPACK_ASSERT(is.good());
    string line;
    while (std::getline(is, line) && is.good()) {
        if (line.find(s) != line.npos) return is;
    }
    throw std::runtime_error("Reached end of file while looking for: " + s);
}

/******************************************************************************************/

std::istream & skip_comments(std::istream &is, string const &comment_start) {
    NUPACK_ASSERT(is.good());
    string line;
    while (is_on_next_line(is, comment_start)) {
        if (!is.good()) throw std::runtime_error("Reached end of file while skipping comments");
        std::getline(is, line);
    }
    return is;
}

/******************************************************************************************/

string peek(std::istream &is){
    string line;
    auto cur = is.tellg();
    std::getline(is, line);
    is.seekg(cur);
    return line;
}

/******************************************************************************************/

std::size_t repeat_char_length(std::string_view s) noexcept {
    std::size_t out = 0;
    auto const sep = sequence_separator();
    for (char const *b=s.data(), *e = b + s.size(); b != e;) {
        if (sep(*b)) {
            ++b;
        } else if (std::isdigit(*b)) {
            out += std::strtoul(b, const_cast<char **>(&b), 10) - 1u;
        } else {
            ++out;
            ++b;
        }
    }
    return out;
}

/******************************************************************************************/

/// Convert dot-parens to pair array
void to_pairs(iseq *v, iseq const n, char const *dp) {
    auto const sep = sequence_separator();
    iseq i = 0, s = n - 1;
    // Number of times to repeat a character, 1 if not given
    auto repeat = [](auto &c) {
        auto d = std::strtoul(++c, const_cast<char **>(&c), 10);
        return d ? d : 1u;
    };
    char const *c = dp;
    while (*c) {
        if (*c == '(') {
            for (auto d = repeat(c); d--; ++i) v[s--] = i; // push on stack
        } else if (*c == ')') {
            auto d = repeat(c);
            if (s + d >= n) NUPACK_ERROR("unmatched ) parenthesis", dp, i);
            for (; d--; ++i) {
                v[i] = v[++s]; // pop off stack
                v[v[i]] = i;
            }
        } else if (*c == '.') {
            for (auto d = repeat(c); d--; ++i) v[i] = i;
        } else if (sep(*c)) {
            ++c;
        } else {
            auto index = c - dp;
            NUPACK_ERROR("bad dot-parens character", dp, index, *c, int(*c));
        }
    };
    if (s != n - 1) NUPACK_ERROR("unmatched ( parenthesis");
    if (i != n) NUPACK_ERROR("dot-parens-plus parsing failed");
}

}

/******************************************************************************************/

Sequence::Sequence(std::string_view letters) {
    // NUPACK_ASSERT(!letters.empty(), "sequence should contain one or more nucleotides");
    auto const n = io::repeat_char_length(letters);
    base_type::resize(n);
    auto const e = io::repeat_char(base_type::begin(), letters);
    if (e != base_type::end())
        NUPACK_ERROR("invalid nucleic acid sequence", letters);
}


Sequence reverse_complement(Sequence seq) noexcept {
    std::reverse(begin_of(seq), end_of(seq));
    for (auto &s : seq) s = complement(s);
    return seq;
}


Sequence reverse_wobble_complement(Sequence seq) noexcept {
    std::reverse(begin_of(seq), end_of(seq));
    for (auto &s : seq) s = wobble_complement(s);
    return seq;
}

/******************************************************************************************/

vec<string> split_sequence_string(string_view s2) {
    vec<string> strs;
    string s(s2);
    boost::trim_if(s, sequence_separator());
    boost::split(strs, s, sequence_separator(), boost::token_compress_on);
    return strs;
}


/******************************************************************************************/

vec<std::array<typename PairList::value_type, 4>> PairList::pseudoknots() const {
    vec<std::array<value_type, 4>> out;
    for_pseudoknots(*this, [&](auto ...is) {out.push_back({static_cast<value_type>(is)...});});
    return out;
}

/******************************************************************************************/

PairList PairList::with_null_bases(small_vec<iseq> const &strand_lengths) const {
    constexpr auto null = std::numeric_limits<iseq>::max();
    data_type out(len(values) + 2 * len(strand_lengths), null);
    // Copy in old values, skip null bases
    auto d = values.begin();
    auto o = out.begin() + 1;
    for (auto l : strand_lengths) {
        std::copy_n(d, l, o);
        d += l;
        o += l + 2;
    }
    // Offset values
    iseq offset = 2 * len(strand_lengths) - 1, end = len(values);
    for (auto l : reversed(strand_lengths)) {
        auto next_end = end - l;
        for (auto &i : out) if (i >= next_end && i < end) i += offset;
        end = next_end;
        offset -= 2;
    }
    // Fix up null bases
    izip(out, [](auto i, auto &j) {if (j == null) j = i;});
    if (!Release) izip(out, [&](auto i, auto &j) {NUPACK_REQUIRE(out[j], ==, i, values, out);});
    return PairList{std::move(out)};
}

/******************************************************************************************/

std::size_t PairList::symmetry() const {
    auto v = values;
    // vector of circular offsets
    izip(v, [s=v.size()](auto i, auto &j) {j = i < j ? j - i : s + j - i;});
    return rotational_symmetry(v);
}

/******************************************************************************************/

}

/******************************************************************************************/

namespace std {
    size_t hash<nupack::Sequence>::operator()(nupack::Sequence const &s) const {
        return nupack::range_hash(nupack::ptr_view(&(s.data()->value), s.size()));
    }

    size_t hash<nupack::Strand>::operator()(nupack::Strand const &s) const {
        return hash<nupack::Sequence>()(s);
    }

    size_t hash<nupack::Complex>::operator()(nupack::Complex const &m) const {
        return nupack::hash_of(m.catenated, m.positions);
    }
}

