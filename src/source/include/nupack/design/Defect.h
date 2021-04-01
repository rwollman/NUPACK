#pragma once
#include "TypeImports.h"

#include "../types/Structure.h"
#include "../types/Matrix.h"

namespace nupack { namespace newdesign {


using defect_vec = vec<std::pair<uint, real>>;

/** @brief simple structure that can represent a per-nucleotide defect at any
  level (complex, tube, multitube). With DesignSequence, can be transformed
  into weights in same order as associated variables. */
struct Defect {
    defect_vec contributions;

    Defect() = default;
    Defect(defect_vec contribs) : contributions(std::move(contribs)) {}
    Defect(vec<real> const &defs, real normalization=1);
    auto total() const {return sum(indirect_view(contributions, [] (auto const &d) {return d.second;}));}
    bool is_valid() {return total() >= 0.0;}

    Defect reduced() const;
    Defect weighted(vec<real> const &weights) const;
    Defect scaled(real) const;

    vec<uint> sample_nucleotides(uint num=1) const;
    // vec<uint> sample_nucleotides(real_mat) const;

    NUPACK_REFLECT(Defect, contributions);
};

vec<real> nucleotide_defects(ProbabilityMatrix const &, Structure const &);
vec<real> nucleotide_defects(Tensor<real, 2> const &, Structure const &);

}}
