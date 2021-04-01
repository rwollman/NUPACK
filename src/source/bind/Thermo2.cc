#include "Thermo.h"
#include <nupack/thermo/CachedModel.h>
#include <nupack/thermo/Adapters.h>
#include <nupack/thermo/ComplexSampler.h>
#include <nupack/Forward.h>
#include <nupack/model/Model.h>

namespace nupack::thermo {

/******************************************************************************************/

void render(Document &doc, Type<ComplexSampler> t) {
    doc.type(t, "thermo.ComplexSampler");
    doc.method(t, "new", rebind::construct<StrandList const &, vec<std::pair<Complex, real>> const &, real>(t));
    doc.method(t, "()", [](ComplexSampler &s, Local &env, CachedModel<PF, Model<>> &mod, usize n) {return s(env, mod, n);});
}

void render_pf(Document &doc) {
    render_lru<3, real32, real64, overflow<real32>>(doc);
    render_engine<PF, 3, 0, 0, 1>(doc, pack<real32, real64, real32>(), as_pack<EnsembleType>());

    render_lru<3, real64, overflow<real32>>(doc);
    render_engine<PF, 3, 0, 1>(doc, pack<real64, real32>(), as_pack<EnsembleType>());

    render_lru<3, real64>(doc);
    render_engine<PF, 3, 0>(doc, pack<real64>(), as_pack<EnsembleType>());

    render_lru<3, overflow<real32>>(doc);
    render_engine<PF, 3, 1>(doc, pack<real32>(), as_pack<EnsembleType>());

    render_lru<3, overflow<real64>>(doc);
    render_engine<PF, 3, 1>(doc, pack<real64>(), as_pack<EnsembleType>());
    doc.render<ComplexSampler>();
}

/******************************************************************************************/

}

