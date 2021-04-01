#pragma once
#include <string>
#include <string_view>

#include "Sequence.h"

namespace nupack {

/******************************************************************************************/

template <bool Subview>
struct ComplexView : ConstIndexable<ComplexView<Subview>>, MemberOrdered {
    small_vec<int> positions; // ending position of each sequence (0 is not included)
    if_t<Subview, View<const_iterator_of<Strand>>, Strand> catenated;
    int offset = 0;
    NUPACK_REFLECT(ComplexView, positions, catenated, offset);

    template <class T, class U>
    ComplexView(T &&t, U &&u, int n=0) : positions(fw<U>(u)), catenated(fw<T>(t)), offset(n) {}

    ComplexView() = default;

    auto const & iter() const {return catenated;}
    auto n_strands() const {return len(positions);}

    auto views() const {
        return vmap<default_vec>(indices(positions), [&](auto i) {return view(catenated, (i ? positions[i-1] : 0), positions[i]);});
    }

    iseq length(iseq i) const {return i ? (positions[i] - positions[i-1]) : positions[i];}

    int first_nick() const {return *(begin_of(positions));}
    int last_nick() const {return *(end_of(positions)-2);}
    auto nicks() const {return view(begin_of(positions), end_of(positions)-1);}
    auto strands() const {auto v = views(); return StrandList(v.begin(), v.end());}
    bool multi() const {return n_strands() > 1;}

    auto slice(iseq b, iseq e) const {
        small_vec<int> positions2{view(positions, b, e)};
        if (b) for (auto &i : positions2) i -= positions[b-1];
        auto seq = view(catenated, b ? positions[b-1] : 0, e ? positions[e-1] : 0);
        int offset = b ? positions[b-1] : 0;
        return ComplexView<true>{seq, std::move(positions2), offset};
    }

    auto strands_included(int i, int j) const {
        auto const b = begin_of(positions), e = end_of(positions);
        return slice(std::upper_bound(b, e, i) - b, std::upper_bound(b, e, j) - b+1);
    }

    friend std::ostream & operator<<(std::ostream &os, ComplexView const &t) {dump_os(os, t.views()); return os;}
};

/******************************************************************************************/

class Complex : public ComplexView<false> {
    Complex(Strand &&s, std::size_t n) : base_type{std::move(s), small_vec<int>(1u, n)} {}
public:
    using base_type = typename Complex::ComplexView;

    Complex() = default;

    template <class V=StrandList, NUPACK_IF(!std::is_base_of_v<Complex, V> && can_construct<Strand, value_type_of<V>>)>
    explicit Complex(V const &v, bool b=true) : base_type{join<Strand>(v), prefixes<small_vec<int>>(false, indirect_view(v, len))} {}

    Complex(StrandList const &v) : Complex(v, true) {}

    explicit Complex(Strand s) : Complex{std::move(s), len(s)} {}

    template <bool Subview> Complex(ComplexView<Subview> const &o)
        : base_type(o.catenated, o.positions) {}

    template <bool Subview> Complex(ComplexView<Subview> &&o) noexcept
        : base_type(std::move(o.catenated), std::move(o.positions)) {}

    Complex duplicated(std::size_t n=2) const {
        Complex out;
        out.catenated = ::nupack::duplicate(catenated, n);
        out.positions = ::nupack::duplicate(positions, n);
        for (auto m : range(1, n)) for (auto i : indices(positions))
            out.positions[i + m * len(positions)] += m * len(catenated);
        return out;
    }

    void rotate_lowest() {
        auto const v = views();
        auto const lowest = lowest_rotation(v);
        if (v != lowest) *this = Complex(lowest);
    }

    auto save_repr() const {return strands();}
    void load_repr(decay<decltype(declval<base_type>().strands())> const &s) {*this = Complex(s);}
};

/******************************************************************************************/

SequenceList complex_to_loop(Complex const &c, int nick);

}

namespace std {
    template <> struct hash<nupack::Complex> {
        size_t operator()(nupack::Complex const &m) const;
    };
}
