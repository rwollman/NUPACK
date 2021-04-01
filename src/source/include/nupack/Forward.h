#pragma once
#include "common/Config.h"
#include "algorithms/Traits.h"

namespace nupack {

/**************************************************************************************/

struct ParameterFile;

struct ParameterInfo;

template <class T>
struct ParameterSet;

template <class T>
struct ParameterData;

template <class T=real>
struct Model;

class System;

struct PairList;

struct Structure;

struct StateBase;

class SequenceSet;

struct JumpSequenceSet;

template <class SS=JumpSequenceSet>
struct JumpLoop;

template <class SS=SequenceSet>
struct StaticLoop;

template <class Base_=StateBase, class Loop_=StaticLoop<>>
struct StaticState;

using State = StaticState<>;

template <class Loop_=JumpLoop<>, class EM=Model<real>>
struct JumpState;

struct Local;

/**************************************************************************************/

NUPACK_NAMESPACE(kmc);
namespace kmc {
    struct EnumeratedObserver;
    struct Stopwatch;
    struct Timer;
    struct FirstCollision;
    template <bool B=true> struct TimeIntegrator;
    template <bool B=true> struct ScaledIntegrator;
    struct HammingObserver;
    struct HittingTimeObserver;
    template <class> struct CovarianceIntegrator;
    template <class> class PairIntegrator;
    template <class> class PairProbabilityBoltzmann;
}

/**************************************************************************************/

NUPACK_NAMESPACE(design);
NUPACK_NAMESPACE(newdesign);
NUPACK_NAMESPACE(concentration);
NUPACK_NAMESPACE(lapack);
NUPACK_NAMESPACE(mma);
NUPACK_NAMESPACE(matlab);
NUPACK_NAMESPACE(memory);
NUPACK_NAMESPACE(simd);

/**************************************************************************************/

NUPACK_NAMESPACE(thermo);
namespace thermo {

NUPACK_NAMESPACE(coax);

struct PF;
struct MFE;

template <class Rig, class Model>
class CachedModel;

}

/**************************************************************************************/


}
