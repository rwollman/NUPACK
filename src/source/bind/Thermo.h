#pragma once
#include "Bind.h"
#include <nupack/thermo/Engine.h>
#include <nupack/thermo/CachedModel.h>
#include <nupack/thermo/ComplexSampler.h>
#include <nupack/types/Structure.h>
#include <nupack/execution/Local.h>

namespace nupack::thermo {

void render(Document &, Type<MemoryLimit>);

/******************************************************************************************/

template <class T, std::size_t N, NUPACK_IF(is_overflow<T>)>
std::pair<rebind::ArrayView, rebind::ArrayView> response(std::type_index, Tensor<T, N> const &m) {
    return {{m.storage.first.data(), rebind::ArrayLayout(m.shape(), m.strides())},
            {m.storage.second.data(), rebind::ArrayLayout(m.shape(), m.strides())}};
}

/******************************************************************************************/

template <class T, std::size_t N, NUPACK_IF(std::is_scalar_v<T>)>
rebind::ArrayView response(rebind::TypeIndex, Tensor<T, N> const &m) {
    return {std::addressof(*m.begin()), rebind::ArrayLayout(m.shape(), m.strides())};
}

/******************************************************************************************/

template <class T, std::size_t N>
void render(Document &doc, Type<Tensor<T, N>> t) {doc.type(t, "thermo.Tensor");}

/******************************************************************************************/

template <class L, NUPACK_IF(is_lru<L>)>
void render(Document &doc, Type<L> t) {
    doc.type(t, "core.LRUCache");
    doc.method(t, "new", rebind::construct<typename L::limit_type::length_type>(t));
    doc.method(t, "keys", [](L const &o) {return as_vec(key_view(o));});
    doc.method(t, "values", [](L const &o) {return as_vec(item_view(o));});
    doc.method(t, "clear", &L::clear);
    doc.method(t, ".limit", &L::limit);
    doc.method(t, ".contents", &L::contents);
}

/******************************************************************************************/

inline void render(Document &doc, Type<PairingAction> t) {
    doc.type(t, "thermo.PairingAction");
    // doc.method(t, "new", rebind::construct<rebind::Callback<bool>>(t));
    doc.method(t, "new", rebind::construct(t));
}

/******************************************************************************************/

template <class T>
void render(Document &doc, Type<XTensor<T>> t) {
    doc.type(t, "thermo.XTensor");
    render_public(doc, t);
}

/******************************************************************************************/

template <class T, NUPACK_IF(traits::is_record<T>)>
void render(Document &doc, Type<T> t) {doc.type(t, "thermo.Record");}

template <class T, NUPACK_IF(traits::is_block<T>)>
void render(Document &doc, Type<T> t) {doc.type(t, "thermo.Block");}

/******************************************************************************************/

template <int N, class Ensemble, class ...Ts>
void render(Document &doc, Type<Cache<N, Ensemble, Ts...>> t) {
    using C = Cache<N, Ensemble, Ts...>;
    doc.render<typename C::base_type>();
    doc.type(t, "thermo.Cache", std::make_tuple(N,  EnsembleType(Ensemble()).index(), overflow_bits<Ts>...));
    // <base_type_of<L>>
    doc.method(t, "new", rebind::construct<std::size_t>(t));
    doc.method(t, "[]", [](C const &l, Complex const &k) {
        auto it = l.find(k);
        if (it == end_of(l)) throw std::out_of_range("Key not found in LRU");
        return it->second;
    });
}

/******************************************************************************************/

template <class M, NUPACK_IF(traits::is_parameter_cache<M>)>
void render(Document &doc, Type<M> t) {
    doc.type(t, "thermo.ParameterCache");
    render_public(doc, t);
}

/******************************************************************************************/

template <class M, NUPACK_IF(traits::is_cached_model<M>)>
void render(Document &doc, Type<M> t, int=0) {
    doc.render<typename M::model_type>();
    doc.type(t, "thermo.CachedModel", std::make_pair(CHAR_BIT * sizeof(value_type_of<M>),
                                                     std::is_same_v<typename M::rig_type, PF>)); // <base_type_of<M>>
    doc.method(t, "new", rebind::construct<typename M::model_type>(t));
    doc.method(t, "reserve", &M::reserve);
    doc.method(t, "set_beta", &M::set_beta);
    doc.method(t, "boltz", &M::template boltz<true>);
    doc.method(t, "capacity", &M::capacity);
    render_public(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<CachedModel<MFE, Model<real32>>> t);
void render(Document &doc, Type<CachedModel<PF, Model<real64>>> t);
void render(Document &doc, Type<CachedModel<PF, Model<real32>>> t);

/******************************************************************************************/

template <class M>
void render(Document &doc, Type<CoaxialRows<M>>) {}

void render(Document &doc, Type<ComplexSampler> t);

/******************************************************************************************/

template <class T, NUPACK_IF(traits::is_message<T>)>
void render(Document &doc, Type<T> t) {
    doc.type(t, "thermo.Message");
    render_public(doc, t);
    doc.method(t, "bits", [](T const &) {return overflow_bits<value_type_of<T>>;});
}

/******************************************************************************************/

template <class D>
std::true_type check_cache_dangle(D, real);

template <class D, int N, class D2, class ...Ts>
std::is_same<D2, D> check_cache_dangle(D, Cache<N, D2, Ts...>);

template <class Rig, int N, bool ...Bs, class ...Types, class ...Dangles>
void render_engine(rebind::Document &doc, pack<Types...> ts, pack<Dangles...> ds) {
    using Caches = rebind::Pack<real, Cache<N, Dangles, oflow<Bs, Types>...> &...>;
    using Models = std::tuple<CachedModel<Rig, Model<Types>> &...>;

    pack<oflow<Bs, Types>...>::for_each([&doc](auto t) {
        NUPACK_UNPACK(doc.render<Message<decltype(*t), Dangles, N>>());
        doc.render<XTensor<decltype(*t)>>();
    });

    using Obs = rebind::Callback<void>;
    using boolCall = rebind::Callback<bool>;

    Caches::for_each([&doc](auto cache) {
        using C = decltype(*cache);
        doc.function("thermo.dynamic_program", [](rebind::Caller call, Local env, Complex const &cx, Models m, C c, Obs o, PairingAction const &a) {
            return dynamic_program<N, Bs...>(env, cx, m, c, std::move(o), a);
        });

        doc.function("thermo.pair_probability", [](rebind::Caller call, Local env, Complex const &cx, Models m, C c, Obs o, PairingAction const &a) {
            return pair_probability<N, Bs...>(env, cx, m, c, std::move(o), a);
        });

        doc.function("thermo.permutations", [](Local env, usize n, Complex const &cx, Models m, C c, Obs o, PairingAction const &a) {
            return permutations<N, Bs...>(env, n, cx.strands(), m, c, std::move(o), a);
        });

        if constexpr(std::is_same_v<Rig, PF>) {
            doc.function("thermo.sample", [](Local env, usize n, usize m, Complex const &cx, Models ms, C c, Obs o, PairingAction const &a) {
                return sample<N, Bs...>(env, n, m, cx, ms, c, std::move(o), a);
            });
        } else {
            doc.function("thermo.subopt", [](Local env, float gap, Complex const &cx, Models m, C c, Obs o, PairingAction const &a, bool print_segments) {
                auto vec = subopt<Outer_Stack, N, Bs...>(env, gap, cx, m, c, std::move(o), a, print_segments);
                return unique_subopt(std::move(vec), cx, std::get<0>(m).energy_model);
            });
            doc.function("thermo.subopt_stream", [](Local env, float gap, Complex const &cx, Models m, C c, boolCall cb, Obs o, PairingAction const &a) {
                auto wrap = [&](auto it) {return cb(*it);};
                subopt_stream<Outer_Stack, N, Bs...>(env, gap, cx, m, wrap, c, std::move(o), a);
            });
        }

        doc.function("thermo.block", [](Local env, Complex const &cx, Models m, C c, Obs o, PairingAction const &a) {
            rebind::Sequence out;
            fork(first_of(m).energy_model.ensemble_type(), [&](auto d) {
                if constexpr(decltype(check_cache_dangle(d, c))::value) {
                    auto b = block<N, Bs...>(env, d, cx, m, c, std::move(o), a);
                    out.emplace_back(b.second);
                    fork(b.first, [&](auto &block) {
                        out.emplace_back(block.names());
                        for_each(block.members(), [&](auto &m) {out.emplace_back(std::move(m));});
                    });
                }
                else throw std::runtime_error("Cache and Model types do not have same dangle setting");
            });
            return out;
        });
    });
}

/******************************************************************************************/

template <int N, class... Ts>
void render_lru(Document &doc) {
    for (auto e : AllEnsembles) fork(ensemble_variant(e), [&](auto d) {
        static_assert(traits::is_cache<Cache<N, decltype(d), Ts...>>, "");
        doc.render<Cache<N, decltype(d), Ts...>>();
    });
}

/******************************************************************************************/

}
