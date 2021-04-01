#include <nupack/types/Structure.h>
#include <nupack/reflect/Serialize.h>

#include <regex>

namespace nupack {

/**
 * @brief Convert a hybrid dpp/dpp-rle structure specification string into a
 *     pure dpp string
 * @details Takes the input string, which is either in dot-parens-plus (dpp)
 *     (e.g. (((+....))) ), dot-parens-plus run-length-encoded (dpp-rle) (e.g.
 *     (3+.4)3 ), or any hybrid of the two (e.g. ((2+...2)3 ), and converts
 *     this string into the pure dpp representation. Throws if the input
 *     string is not a hybrid (including pure dpp or dpp-rle).
 *
 * @param s the hybrid string to convert to pure dpp representation
 * @return a pure dpp representation of the structure implied by s.
 */
string Structure::parse_struc(string_view s0) {
    string s(s0); // fix sometime...
    /* Return simple dpp structure */
    static std::regex pure_dpp("[().+]*");
    if (std::regex_match(s, pure_dpp)) return s;

    /* If not a simple dpp structure, make sure it's a valid hybrid */
    static std::regex rle_dpp_whole("([().][0-9]*|[+])+");
    if (!std::regex_match(s, rle_dpp_whole)) throw std::invalid_argument(s + " is not in dpp or rle_dpp format");

    static std::regex rle_dpp_component("([().])([0-9]*)|([+])");
    auto beg = std::sregex_iterator(s.begin(), s.end(), rle_dpp_component);
    auto end = std::sregex_iterator();

    std::stringstream ss;
    for (auto i = beg; i != end; ++i) {
        auto m = indirect_view(view(i->begin(), i->end()), [](auto const &j) {return j.str();});
        int repeats = m[2].empty() ? 1 : std::stoi(m[2]);
        string cur = m[3].empty() ? m[1] : m[3];
        for (auto i : range(repeats)) ss << cur;
    }

    return ss.str();
}

/******************************************************************************************/

/**
 * @brief Return a minimal length run-length enconding of the dot-parens-plus structure
 * @details For single character runs, e.g. the unpaired base in (((+.))), the
 *     number "1" will not be included in the output string as it is
 *     redundant.
 * @return a minimal length dpp-rle string representation of the Structure.
 */
string Structure::dp_rle() const {
    auto s = dp();
    static std::regex dp_run("\\(+|\\.+|\\)+|\\+");

    std::stringstream ss;
    auto end = std::sregex_iterator();
    for (auto i = std::sregex_iterator(s.begin(), s.end(), dp_run); i != end; ++i) {
        auto run = i->str();
        auto n = len(run);
        if (n == 1) ss << run;
        else ss << run[0] << std::to_string(n);
    }

    return ss.str();
}

/******************************************************************************************/

std::size_t Structure::symmetry() const {
    auto n = nicks;
    std::adjacent_difference(n.begin(), n.end(), n.begin());
    auto sym1 = rotational_symmetry(n);
    return sym1 == 1 ? 1 : std::gcd(sym1, PairList::symmetry());
}

/******************************************************************************************/

void PairList::rotate(std::ptrdiff_t s) {
    izip(values, [&](auto i, auto &j) {j = (j - s) % values.size();});
    std::rotate(values.begin(), values.begin() + s % values.size(), values.end());
}

void Structure::rotate(std::ptrdiff_t s) {
    if (!s) return;
    s %= values.size();
    PairList::rotate(nicks[s-1]);
    std::adjacent_difference(nicks.begin(), nicks.end(), nicks.begin()); // strand lengths
    std::rotate(nicks.begin(), nicks.begin() + s, nicks.end());
    std::partial_sum(nicks.begin(), nicks.end(), nicks.begin()); // now redo
}

void PairList::append(PairList const &p) {
    auto const n = len(values);
    values.reserve(n + len(p));
    for (auto i : p) values.emplace_back(n + i);
}

void PairList::throw_if_invalid() const {
    NUPACK_REQUIRE(len(values), >, 0, "Empty pair list");
    izip(values, [&](std::size_t i, std::size_t j) {
        NUPACK_REQUIRE(j, <, len(values), "Pair index too large", values);
        NUPACK_REQUIRE(i, ==, values[j], "Mismatched base pair", values, i, j);
    });
}

bool PairList::is_connected(small_vec<uint> const &nicks) const {
    NUPACK_ASSERT(!values.empty(), "empty pairs");
    NUPACK_ASSERT(!nicks.empty(), "empty nicks");

    small_vec<bool> visited(len(values));
    std::vector<std::pair<uint, uint>> stack;

    stack.emplace_back(0, 0);

    while (!stack.empty()) {
        auto [i, s] = stack.back();
        stack.pop_back();
        visited[i] = true;

        if (auto j = values[i]; !visited[j]) {
            stack.emplace_back(j, std::upper_bound(nicks.begin(), nicks.end(), j) - nicks.begin());
        }

        if (i+1 < nicks[s] && !visited[i+1]) {
            stack.emplace_back(i+1, s);
        }

        if (i > (s ? nicks[s-1] : 0) && !visited[i-1]) {
            stack.emplace_back(i-1, s);
        }
    }

    return all_of(visited);
}

/******************************************************************************************/

}

namespace std {

size_t hash<nupack::PairList>::operator()(nupack::PairList const &p) const {return nupack::hash_of(p.values);}

size_t hash<nupack::Structure>::operator()(nupack::Structure const &p) const {return nupack::hash_of(p.nicks, p.values);}

}
