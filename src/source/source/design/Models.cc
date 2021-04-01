#include <nupack/design/Models.h>

namespace nupack { namespace newdesign {


void ThermoEnviron::initialize_cache(std::size_t ram) {
    fork(std::get<0>(models).energy_model.ensemble_type(), [&](auto x) {
        using cache_type = DesignCache<decltype(x)>;
        cache = cache_type(ram);
    });
}


ThermoEnviron & ModelMap::get(Model<real> const &key) const {
    // std::lock_guard<std::mutex> lk(mut.mut)
    return mod_map.try_emplace(key, key).first->second;
}


void ThermoEnviron::add_pfunc(::nupack::Complex const &s, real log_pfunc) {
    std::unique_lock lock(mut.mut);
    log_pfuncs.try_emplace(s, log_pfunc);
}


Variant<bool, real> ThermoEnviron::get_pfunc(::nupack::Complex const &s) const {
    std::shared_lock lock(mut.mut);
    auto it = log_pfuncs.find(s);
    if (it == end_of(log_pfuncs)) return false;
    return it->second;
}


void ThermoEnviron::clear_cache() {
    fork(cache, [&](auto &x) {x.clear();});
}


void ModelMap::clear_caches() {
    for (auto &i : item_view(mod_map)) i.clear_cache();
}


void ModelMap::create_caches(std::size_t ram) {
    auto num_models = len(mod_map);
    auto ram_per_model = ram / num_models;
    for (auto &i : item_view(mod_map)) i.initialize_cache(ram_per_model);
}

std::size_t ModelMap::ram() const {
    std::size_t total {0};
    for (auto &i: item_view(mod_map)) {
        total += fork(i.cache, [&](auto const &x) {return x.limit.length;});
    }
    return total;
}

}}
