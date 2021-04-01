#pragma once
#include "Bind.h"
#include <nupack/model/Model.h>

namespace rebind {

template <> struct ImplicitConversions<nupack::Model<float>> {using types = Pack<nupack::Model<double>>;};
template <> struct ImplicitConversions<nupack::Model<double>> {using types = Pack<nupack::Model<float>>;};

}

/******************************************************************************************/

namespace nupack {

/******************************************************************************************/

void render(Document &doc, Type<Model<real32>>);
void render(Document &doc, Type<Model<real64>>);
void render(Document &doc, Type<EnsembleType>);
void render(Document &doc, Type<ParameterFile>);
void render(Document &doc, Type<ParameterData<real>>);
void render(Document &doc, Type<ParameterData<float>>);
void render(Document &doc, Type<ParameterInfo>);
void render(Document &doc, Type<ParameterSet<real32>>);
void render(Document &doc, Type<ParameterSet<real64>>);
void render(Document &doc, Type<ModelConditions>);

/******************************************************************************************/

template <class M>
auto loop_energy(M const &model, Complex const &c, int nick) {
    auto const v = complex_to_loop(c, nick);
    return model.pairable.check_loop(v) ? model.loop_energy(v, nick) : *inf;
}

/******************************************************************************************/

template <class T>
void render(Document &doc, Type<Model<T>> t, int=0) {
    using M = Model<T>;
    doc.type(t, "model.Model", CHAR_BIT * sizeof(T));
    render_public(doc, t);
    render_comparisons(doc, t);

    doc.method(t, "new", [](Ensemble e, ParameterFile const &p, ModelConditions const &cs, uint gu) {
        auto mod = Model<T>(e, p, cs, gu == 2 ? std::nullopt :
            std::make_optional(gu ? WobblePairing::on : WobblePairing:: off));
        set_sequence_type_weak(mod.parameters.material == "RNA");
        return mod;
    });
    doc.method(t, "join_penalty",         &M::join_penalty);
    doc.method(t, "multi_init",           &M::multi_init);
    doc.method(t, "multi_base",           &M::multi_base);
    doc.method(t, "multi_pair",           &M::multi_pair);
    doc.method(t, "interior_size_energy", &M::interior_size_energy);
    doc.method(t, "interior_asymmetry",   &M::interior_asymmetry);
    doc.method(t, "interior_mismatch",    &M::interior_mismatch);
    // doc.method(t, "dangle5",              &M::dangle5);
    // doc.method(t, "dangle3",              &M::dangle3);
    doc.method(t, "boltz",                &M::boltz);
    doc.method(t, "hairpin_energy",       &M::template hairpin_energy<Sequence>);
    doc.method(t, "loop_energy",          &loop_energy<M>);
    doc.method(t, "stack_energies",       &loop_stacking_energies<M>);
    doc.method(t, "multi_energy",         &M::template multi_energy<SequenceList>);
    doc.method(t, "exterior_energy",      &M::template exterior_energy<SequenceList>);
    doc.method(t, "interior_energy",      &M::template interior_energy<Sequence, Sequence>);
    doc.method(t, "coaxial_stack_energy", &M::coaxial_stack_energy);
    render_json(doc, t);
}

/******************************************************************************************/

template <class T>
void render(Document &doc, Type<ParameterData<T>> t, int=0) {
    doc.type(t, "model.ParameterData");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "new", [](ParameterFile const &file, string const &kind) {return ParameterData<T>(file.open().at(kind));});
    // doc.method(t, "array", [](ParameterData<T> const &p) {
    //     NUPACK_ASSERT(p.array, "empty parameters");
    //     return arma::Col<T>(p.array.get(), ParameterData<T>::size);
    // });
    render_json(doc, t);
    render_public(doc, t);
}

/******************************************************************************************/

template <class T>
void render(Document &doc, Type<ParameterSet<T>> t, int=0) {
    doc.render<ParameterData<T>>();
    doc.type(t, "model.ParameterSet");
    doc.method(t, "new", rebind::construct<ParameterInfo>(t));
    render_json(doc, t);
    render_public(doc, t);
}

/******************************************************************************************/

}
