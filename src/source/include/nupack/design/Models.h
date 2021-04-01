#pragma once

#include "TypeImports.h"
#include "../thermo/CachedModel.h"
#include "../model/Model.h"
#include <shared_mutex>

namespace nupack::newdesign {

template <class D>
using DesignCache = thermo::Cache<3, D, real32, real64, overflow<real32>, overflow<real64>>;

using DesignCacheVariant = Variant<DesignCache<NoStacking>, DesignCache<MinDangles>, DesignCache<AllDangles>, DesignCache<Stacking>>;

/**
 * @brief interface adapter that "doubles" the length of a set of models to
 * support non-overflow and overflow versions of the thermo code and allow
 * fallback while still only using 2 underlying models
 *
 * @tparam T type of paired 32- and 64-bit models
 * @param t the paired models
 * @return tuple of references to the underlying models in the correct doubled
 * order
 */
template <class T>
auto double_models(T &&t) {
    return std::tie(std::get<0>(t), std::get<1>(t), std::get<0>(t), std::get<1>(t));
};


using ModelsTuple = std::tuple<CachedModel<PF, Model<real32>>, CachedModel<PF, Model<real64>>>;

struct CopyableMutex {
    CopyableMutex() = default;
    CopyableMutex(CopyableMutex const &) : CopyableMutex() {};
    CopyableMutex(CopyableMutex &&) : CopyableMutex() {};
    CopyableMutex & operator=(CopyableMutex const &) { return *this; };
    CopyableMutex & operator=(CopyableMutex &&) { return *this; };

    std::shared_mutex mutable mut;
};

struct ThermoEnviron {
    ModelsTuple models;
    DesignCacheVariant cache;
    std::map<::nupack::Complex, real> log_pfuncs;
    CopyableMutex mut;

    ThermoEnviron() = default;
    /* create cache of correct type immediately when environment is made */
    // ThermoEnviron(ModelsTuple mods) : models(std::move(mods)) {initialize_cache(0);};
    ThermoEnviron(Model<real> const &model) : models(model, model) {initialize_cache(0);};

    void initialize_cache(std::size_t ram);
    auto doubled() const { return double_models(models); }
    void add_pfunc(::nupack::Complex const &s, real log_pfunc);
    Variant<bool, real> get_pfunc(::nupack::Complex const &s) const;
    void clear_cache();

    NUPACK_REFLECT(ThermoEnviron, models, cache, log_pfuncs);
};


/**
 * @brief Maintains cache of all models needed during design to avoid constantly
 * recreating the model. Furthermore, keeps 32-bit and 64-bit versions of the
 * same model together to support model fallback operations in a relatively
 * seamless manner when calling thermo code.
 *
 */
class ModelMap : MemberOrdered {
    mutable std::map<Model<real>, ThermoEnviron> mod_map;

public:
    NUPACK_REFLECT(ModelMap, mod_map);
    ModelMap() = default;

    /**
     * @brief Create new model if key has not been seen before and return
     * pointer to model.
     *
     * @param key specification of the requested model
     * @return reference to the requested pair of models
     */
    ThermoEnviron & get(Model<real> const & key) const;
    auto const &cached_models(Model<real> const & key) const {return get(key).models;}

    void create_caches(std::size_t ram);
    void clear_caches();
    auto size() const {return len(mod_map);}
    std::size_t ram() const;
};



}
