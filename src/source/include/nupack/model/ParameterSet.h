#pragma once
#include "../standard/Array.h" // for multi_array
#include "../common/Constants.h" // for temperature
#include "../types/IO.h" // for goto_line_after
#include "../types/Sequence.h"
#include "../types/Matrix.h"
#include "../standard/Ptr.h"
#include "../reflect/Serialize.h"
#include "ParameterStorage.h"


namespace nupack {

/******************************************************************************************/

struct ParameterFile : MemberOrdered {

    string path;
    NUPACK_REFLECT(ParameterFile, path);

    /* Take a parameter file name and load in the parameters
    Priority ranking:
        1) the specified path
        2) "NUPACKHOME" environment variable
        3) then the NUPACK parameters folder (set by CMake)
    */
    ParameterFile(string name="rna");

    json open() const;
};

NUPACK_DEFINE_TYPE(is_parameter_file, ParameterFile);

/******************************************************************************************/

template <class T>
struct ParameterData : TotallyOrdered {
    using value_type = T;
    static constexpr std::size_t size = decltype(join_penalty)::end;

    std::shared_ptr<T> array; // should be T[] but seems screwed up to compile on my clang
    static auto allocate() {
        return std::shared_ptr<T>(new T[size], [](T *t) {delete[] t;});
    }

    NUPACK_REFLECT(ParameterData, array);

    T const *begin() const {return array.get();}
    T const *end() const {return array.get() + size;}

    T *begin() {return array.get();}
    T *end() {return array.get() + size;}

    bool operator==(ParameterData const &p) const {
        if (array == p.array || !array || !p.array) return array == p.array;
        return std::equal(begin(), end(), p.begin(), p.end());
    }

    bool operator<(ParameterData const &p) const {
        if (array == p.array || !array || !p.array) return array < p.array;
        return std::lexicographical_compare(begin(), end(), p.begin(), p.end());
    }

    void load_repr(json const &j);
    json save_repr() const;

    ParameterData() = default;
    explicit ParameterData(std::shared_ptr<T> ptr) : array(std::move(ptr)) {}
    explicit ParameterData(json const &j) {load_repr(j);}

    template <class U>
    ParameterData(ParameterData<U> const &p) : array(allocate()) {
        std::copy(p.begin(), p.begin() + size, begin());
    }

    template <class I>
    auto span(I) const {return view(begin() + I::begin, begin() + I::end);}

    template <class I>
    auto span(I) {return view(begin() + I::begin, begin() + I::end);}

    template <class I, class ...Is>
    T const &operator()(I, Is ...is) const {
        auto idx = I::index(is...);
        NUPACK_DASSERT(array, "array is not allocated");
        return begin()[idx];
    }

    void add_loop_bias(T t) {
        for_each(std::tie(stack, bulge_size, interior_size, hairpin_size,
            interior_1_1, interior_1_2, interior_2_2, join_penalty, multi_init),
                [&](auto key) {for (auto &x : span(key)) x += t;});
    }
};

/******************************************************************************************/

struct ParameterInfo {
    ParameterFile file;
    string kind = "dG";
    real loop_bias = 0;
    real temperature = DefaultTemperature;

    NUPACK_REFLECT(ParameterInfo, temperature, loop_bias, kind, file);
    using is_member_ordered = True;
};

/******************************************************************************************/

template <class T>
struct ParameterSet : TotallyOrdered {
    using value_type = T;
    ParameterInfo info;
    ParameterData<T> data;
    std::string material;
    bool default_wobble_pairing;
    NUPACK_REFLECT(ParameterSet, info, data, material, default_wobble_pairing);

    template <class U>
    ParameterSet(ParameterSet<U> const &o)
        : info{o.info}, data(o.data), material(o.material), default_wobble_pairing(o.default_wobble_pairing) {
        static_assert(!is_same<T, U>, "Should use normal copy constructor");
    }

    ParameterSet() = default;

    explicit ParameterSet(ParameterInfo i);

    bool operator==(ParameterSet const &p) const {return info == p.info;}
    bool operator<(ParameterSet const &p) const {return info < p.info;}

    friend std::ostream & operator<<(std::ostream &os, ParameterSet const &p) {
        return os << "ParameterSet(" << p.info.file.path << ", " << p.info.kind << ", " << p.info.temperature << " K)";
    }

    void load_repr(ParameterInfo const &i) {*this = ParameterSet(i);}
    ParameterInfo save_repr() const {return info;}
};

/******************************************************************************************/

template <class T>
ParameterSet<T>::ParameterSet(ParameterInfo i) : info{std::move(i)} {
    json const j = info.file.open();
    data = ParameterData<T>(j.at("dG"));
    j.at("default_wobble_pairing").get_to(default_wobble_pairing);
    j.at("material").get_to(material);

    T const t = info.temperature;
    if (t != DefaultTemperature) {
        T kg = t / DefaultTemperature, kh = 1 - kg;
        zip(data, ParameterData<T>(j.at("dH")), [kg, kh](T &g, T const &h) {
            g = kg * g + kh * h;
        });
    }

    data.begin()[decltype(join_penalty)::begin] -= std::log(water_molarity(t)) * Kb * t;
    data.add_loop_bias(info.loop_bias);
}

/******************************************************************************************/

namespace base_index {

/******************************************************************************************/

// String of bases to an array of indices (known length)
template <std::size_t N, class S>
std::array<BaseIndex, N> nd_sequence_index(S const &key) {
    std::array<BaseIndex, N> out;
    NUPACK_REQUIRE(N, ==, key.size(), "Incorrect number of nucleotides");
    zip(out, key, [](auto &o, auto c) {o = Base::lookup(c);});
    return out;
}

/******************************************************************************************/

template <class V, class T>
void load_array(V &&v, T, json const &x) {
    for (auto const &[key, value] : x.items())
        value.get_to(v[T::array_index(nd_sequence_index<T::ndim>(key))]);
}

template <class ...Is>
std::string nd_string(Is const ...is) {
    std::string out;
    (out.push_back(Base::names[is]), ...);
    return out;
}

/******************************************************************************************/

template <class T, class V, class ...Is>
void save_to_array(json &j, T, V const &v, Is const ...is) {
    if constexpr(sizeof...(Is) == T::ndim) {
        auto value = v[T::index(is...)];
        if (value != 0) j[nd_string(is...)] = value;
    } else {
        for (auto i : range(4)) save_to_array(j, T(), v, i, is...);
    }
}

template <class V, class T>
json save_array(V const &v, T t) {
    json out = json::object();
    save_to_array(out, t, v);
    return out;
}

/******************************************************************************************/

template <class P, class T>
void simple_load(P &p, T, json const &j) {
    if constexpr(T::ndim == 0) {
        j.get_to(*p.span(T()).begin());
    } else {
        auto s = p.span(T());
        NUPACK_REQUIRE(j.size(), ==, s.size());
        std::copy(j.begin(), j.end(), s.begin());
    }
}

template <class P, class T>
json simple_save(P const &p, T) {
    if constexpr(T::ndim == 0) {
        return p(T());
    } else {
        auto out = json::array();
        for (auto const &x : p.span(T())) out.emplace_back(x);
        return out;
    }
}

}

template <class T>
void ParameterData<T>::load_repr(json const &j) {
    if (!array) array = allocate();
    fill(*this, T(0));
    base_index::simple_load(*this, log_loop_penalty,  j.at("log_loop_penalty"));
    base_index::simple_load(*this, hairpin_size,  j.at("hairpin_size"));
    base_index::simple_load(*this, bulge_size,  j.at("bulge_size"));
    base_index::simple_load(*this, multi_init,  j.at("multiloop_init"));
    base_index::simple_load(*this, multi_pair,  j.at("multiloop_pair"));
    base_index::simple_load(*this, multi_base,  j.at("multiloop_base"));
    base_index::simple_load(*this, join_penalty,  j.at("join_penalty"));
    base_index::simple_load(*this, interior_size,  j.at("interior_size"));
    base_index::simple_load(*this, ninio,  j.at("asymmetry_ninio"));

    base_index::load_array(begin(), stack, j.at("stack"));
    base_index::load_array(begin(), coaxial_stack, j.at("coaxial_stack"));
    base_index::load_array(begin(), hairpin_tri, j.at("hairpin_triloop"));
    base_index::load_array(begin(), hairpin_tetra, j.at("hairpin_tetraloop"));
    base_index::load_array(begin(), hairpin_mismatch, j.at("hairpin_mismatch"));
    base_index::load_array(begin(), interior_mismatch, j.at("interior_mismatch"));
    base_index::load_array(begin(), terminal_mismatch, j.at("terminal_mismatch"));
    base_index::load_array(begin(), dangle5, j.at("dangle_5"));
    base_index::load_array(begin(), dangle3, j.at("dangle_3"));
    base_index::load_array(begin(), interior_1_1, j.at("interior_1_1"));
    base_index::load_array(begin(), interior_1_2, j.at("interior_1_2"));
    base_index::load_array(begin(), interior_2_2, j.at("interior_2_2"));
    base_index::load_array(begin(), terminal_penalty, j.at("terminal_penalty"));
}

template <class T>
json ParameterData<T>::save_repr() const {
    json j;
    if (!begin()) return j;
    j["log_loop_penalty"] = base_index::simple_save(*this, log_loop_penalty);
    j["hairpin_size"] = base_index::simple_save(*this, hairpin_size);
    j["bulge_size"] = base_index::simple_save(*this, bulge_size);
    j["multiloop_init"] = base_index::simple_save(*this, multi_init);
    j["multiloop_pair"] = base_index::simple_save(*this, multi_pair);
    j["multiloop_base"] = base_index::simple_save(*this, multi_base);
    j["join_penalty"] = base_index::simple_save(*this, join_penalty);
    j["interior_size"] = base_index::simple_save(*this, interior_size);
    j["asymmetry_ninio"] = base_index::simple_save(*this, ninio);

    j["stack"] = base_index::save_array(begin(), stack);
    j["coaxial_stack"] = base_index::save_array(begin(), coaxial_stack);
    j["hairpin_triloop"] = base_index::save_array(begin(), hairpin_tri);
    j["hairpin_tetraloop"] = base_index::save_array(begin(), hairpin_tetra);
    j["hairpin_mismatch"] = base_index::save_array(begin(), hairpin_mismatch);
    j["interior_mismatch"] = base_index::save_array(begin(), interior_mismatch);
    j["terminal_mismatch"] = base_index::save_array(begin(), terminal_mismatch);
    j["dangle_5"] = base_index::save_array(begin(), dangle5);
    j["dangle_3"] = base_index::save_array(begin(), dangle3);
    j["interior_1_1"] = base_index::save_array(begin(), interior_1_1);
    j["interior_1_2"] = base_index::save_array(begin(), interior_1_2);
    j["interior_2_2"] = base_index::save_array(begin(), interior_2_2);
    j["terminal_penalty"] = base_index::save_array(begin(), terminal_penalty);
    return j;
}

/******************************************************************************************/

}
