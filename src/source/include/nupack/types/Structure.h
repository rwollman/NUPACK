#pragma once
#include "PairList.h"
#include "Sequence.h"

namespace nupack {

/******************************************************************************************/

/** @brief A secondary structure with nicks embedded */
struct Structure : PairList {
    using base_type = PairList;

    Nicks nicks;
    NUPACK_REFLECT_BASE(Structure, PairList, nicks);

    Structure() = default;

    // from DPP
    Structure(string_view s, bool) :
            PairList(s),
            nicks(prefixes<Nicks>(false, indirect_view(split_sequence_string(s), len))) {}

    // from DPP or rle DPP
    explicit Structure(string_view s) : Structure(parse_struc(s), bool()) {}

    Structure(PairList pairs, Nicks nicks) : PairList(std::move(pairs)), nicks(std::move(nicks)) {}

    string dp() const {return PairList::dp(nicks);}
    string dp_rle() const;

    bool is_connected() const {return PairList::is_connected(nicks);}
    bool valid() const {return len(*this) > 0;}

    // Rotational symmetry of this structure
    std::size_t symmetry() const;

    std::size_t n_strands() const {return nicks.size();}

    std::size_t strand_length(std::size_t i) const {return i ? nicks[i] - nicks[i-1] : nicks[0];}

    // Rotate by given number of strands
    void rotate(std::ptrdiff_t i);

    static string parse_struc(string_view s);

    auto save_repr() const {return dp_rle();}
    void load_repr(string s) {*this = Structure(s);}

    friend std::ostream & operator<<(std::ostream &os, Structure const &s) {
        return os << "Structure(\"" << s.dp_rle() << "\")";
    }
};

/******************************************************************************************/

inline Structure operator"" _dp(char const *s) {return Structure(std::string_view(s));}
inline Structure operator"" _dp(char const *s, std::size_t n) {return Structure(std::string_view(s, n));}


/******************************************************************************************/

}


namespace std {
    template <>
    struct hash<nupack::Structure> {
        size_t operator()(nupack::Structure const &p) const;
    };
}
