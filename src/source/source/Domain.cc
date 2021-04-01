#include <nupack/types/Domain.h>
#include <nupack/reflect/Serialize.h>

/******************************************************************************************/

namespace std {
    size_t hash<nupack::NamedComplex>::operator()(nupack::NamedComplex const &c) const {
        return nupack::hash_of(c.lowest_names());
    }
    size_t hash<nupack::TargetComplex>::operator()(nupack::TargetComplex const &c) const {
        return nupack::hash_of(lowest_rotation(c.strands));
    }
}

/******************************************************************************************/

namespace nupack {

TargetComplex::TargetComplex(vec<TargetStrand> v, string name, Structure s, real bonus)
    : strands(std::move(v)), structure(std::move(s)), name(std::move(name)), bonus(bonus) {

    if (!structure.size()) return;

    NUPACK_REQUIRE(len(strands), ==, structure.n_strands(), "TargetComplex(): incorrect number of strands in structure", strands, structure);

    izip(strands, [&](auto i, auto const &s) {
        NUPACK_REQUIRE(s.size(), ==, structure.strand_length(i), "TargetComplex(): incorrect length of strand", i, strands, structure);
    });

    NUPACK_REQUIRE(sum(strands, len), ==, structure.size(), "Sequence and structure sizes do not agree", strands, structure);
}

/******************************************************************************************/

bool TargetComplex::operator==(TargetComplex const &c) const {
    return lowest_rotation(strands) == lowest_rotation(c.strands);
}

/******************************************************************************************/

bool TargetComplex::operator<(TargetComplex const &c) const {
    return lowest_rotation(strands) < lowest_rotation(c.strands);
}

/******************************************************************************************/

Domain Domain::reverse_complement(bool wobble) const {
    auto s = name;
    if (s.empty() || s.back() != '*') s.push_back('*');
    else s.pop_back();
    auto const &seq = static_cast<Sequence const &>(*this);
    return {!complement.empty() ? complement :
        (wobble ? reverse_wobble_complement(seq) : ::nupack::reverse_complement(seq)), seq,
        std::move(s)};
}

/******************************************************************************************/

vec<NamedStrand> NamedComplex::strands() const {
    vec<NamedStrand> out;
    out.reserve(strand_names.size());
    zip(strand_names, Complex::strands(), complements,
        [&](auto const &n, Strand s, Strand c) {out.emplace_back(std::move(s), c, n);});
    return out;
}

/******************************************************************************************/

NamedStrand NamedStrand::reverse_complement(bool wobble) const {
    auto s = name;
    if (s.empty() || s.back() != '*') s.push_back('*');
    else s.pop_back();
    auto const &seq = static_cast<Sequence const &>(*this);
    return {!complement.empty() ? static_cast<Sequence const &>(complement) :
        (wobble ? reverse_wobble_complement(seq) : ::nupack::reverse_complement(seq)), seq,
        std::move(s)};
}

/******************************************************************************************/

vec<std::string_view> NamedComplex::lowest_names() const {
    auto i = lowest_rotational_order(views());
    vec<std::string_view> n(strand_names.begin(), strand_names.end());
    std::rotate(n.begin(), n.begin() + i, n.end());
    return n;
}

/******************************************************************************************/

std::size_t NamedComplex::symmetry() const {
    return std::gcd(rotational_symmetry(strand_names), rotational_symmetry(views()));
}

/******************************************************************************************/

bool NamedComplex::operator<(NamedComplex const &c) const {
    return std::make_tuple(lowest_rotation(  views()),   lowest_names())
         < std::make_tuple(lowest_rotation(c.views()), c.lowest_names());
}

bool NamedComplex::operator==(NamedComplex const &c) const {
    return std::make_tuple(lowest_rotation(  views()),   lowest_names())
        == std::make_tuple(lowest_rotation(c.views()), c.lowest_names());
}

/******************************************************************************************/

}
