#pragma once
#include "BasicBlock.h"
#include "CoaxialBlock.h"
#include "../types/LRU.h"
#include "../standard/Variant.h"

namespace nupack::thermo {

/**************************************************************************************/

/// Type of record emitted by a given Block type
template <class Block>
using BlockRecord = decltype(declval<Block>().write(span(0, 0), span(0, 0), true));

/// Type of record emitted by a Block of the given data type, dangles, and complexity
template <class T, class Ensemble, int N>
using RecordType = BlockRecord<BlockMatrix<T, Ensemble, N>>;

/**************************************************************************************/

/// Convenience wrapper for LRU<...> cache class for easier use with variants
template <int N, class Ensemble, class ...Ts>
struct Cache : LRU<Complex, MaybeVariant<RecordType<Ts, Ensemble, N>...>, MemoryLimit> {
    using base_type = typename Cache::LRU;
    using is_lru = False;
    Cache(std::size_t n=0) : base_type(n) {}
};

NUPACK_DEFINE_VARIADIC(is_cache, Cache, int, class, class);

/**************************************************************************************/

}
