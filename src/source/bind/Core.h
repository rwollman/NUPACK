#pragma once
#include "Bind.h"
#include <nupack/types/Sequence.h>
#include <nupack/types/Domain.h>
#include <nupack/math/Sparse.h>

namespace nupack {

/******************************************************************************/

template <class T>
void render(Document &doc, Type<SparsePairs<T>> t) {
    doc.type(t, "core.SparsePairs");
    doc.method(t, "new", &sparse_pair_matrix<T>);
    render_public(doc, t);
}

/******************************************************************************/

inline auto response(std::type_index, Complex const &v) {return v.strands();}
inline char response(std::type_index, Base const &b) {return b.safe_letter();}
inline string response(std::type_index, Sequence const &s) {return s.str();}

std::optional<Base> request(Type<Base>, rebind::Variable const &r, rebind::Dispatch &msg);
std::optional<Complex> request(Type<Complex> t, rebind::Variable const &v, rebind::Dispatch &msg);
std::optional<PairList> request(Type<PairList>, rebind::Variable const &r, rebind::Dispatch &msg);
std::optional<Sequence> request(Type<Sequence>, rebind::Variable const &r, rebind::Dispatch &msg);
std::optional<Strand> request(Type<Strand>, rebind::Variable const &r, rebind::Dispatch &msg);

void render(Document &doc, Type<AlwaysFalse>);
void render(Document &doc, Type<AlwaysTrue>);
void render(Document &doc, Type<Base> t);
void render(Document &doc, Type<Complex> t);
void render(Document &doc, Type<Domain> t);
void render(Document &doc, Type<Local> t);
void render(Document &doc, Type<NamedComplex> t);
void render(Document &doc, Type<NamedStrand> t);
void render(Document &doc, Type<PairList>);
void render(Document &doc, Type<Pairable> t);
void render(Document &doc, Type<Sequence> t);
void render(Document &doc, Type<Strand> t);
void render(Document &doc, Type<Structure> t);
void render(Document &doc, Type<TargetComplex> t);
void render(Document &doc, Type<TargetStrand> t);

void set_sequence_type_weak(bool rna) noexcept;

}
