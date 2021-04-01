#pragma once

#include "TypeImports.h"
#include "Logging.h"

namespace nupack { namespace newdesign {

struct EngineObserver;
struct ThermoEnviron;

using models_type = std::tuple<CachedModel<PF, Model<real32>> const &, CachedModel<PF, Model<real64>> const &,
                               CachedModel<PF, Model<real32>> const &, CachedModel<PF, Model<real64>> const &>;
using SplitPoint = std::pair<uint, uint>;

real partition_function(Local const &env, ::nupack::Complex const &, ThermoEnviron &, EngineObserver &obs=NullEngineObserver);
std::pair<Tensor<real, 2>, real> pair_probability(Local const &env, ::nupack::Complex const &, ThermoEnviron &, EngineObserver &obs=NullEngineObserver);
real partition_function(Local const &env, ::nupack::Complex const &, models_type const &, EngineObserver &obs=NullEngineObserver);
std::pair<Tensor<real, 2>, real> pair_probability(Local const &env, ::nupack::Complex const &, models_type const &, EngineObserver &obs=NullEngineObserver);
std::pair<Tensor<real, 2>, real> pair_probability(Local const &env, ::nupack::Complex const &, models_type const &, vec<SplitPoint> const &fixed_pairs, real bonus, EngineObserver &obs=NullEngineObserver);

}}
