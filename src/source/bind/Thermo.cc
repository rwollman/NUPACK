#include "Thermo.h"
#include <nupack/thermo/CachedModel.h>
#include <nupack/thermo/Adapters.h>
#include <nupack/Forward.h>
#include <nupack/model/Model.h>


namespace nupack::thermo {

/******************************************************************************************/

template <class T, class Base>
rebind::Variable response(rebind::TypeIndex t, Block<T, Base> const &b) {
    return t.equals<rebind::Dictionary>() ? to_dictionary(b) : rebind::Variable();
}

void render(Document &doc, Type<CachedModel<MFE, Model<real32>>> t) {render(doc, t, 0);}
void render(Document &doc, Type<CachedModel<PF,  Model<real64>>> t) {render(doc, t, 0);}
void render(Document &doc, Type<CachedModel<PF,  Model<real32>>> t) {render(doc, t, 0);}

/******************************************************************************************/

void render(Document &doc, Type<MemoryLimit> t) {
    doc.type(t, "core.MemoryLimit");
    render_public(doc, t);
}

/******************************************************************************************/

void render_pf(Document &doc);

void render_mfe(Document &doc) {
    doc.render<CachedModel<MFE, Model<real32>>>();
    doc.render<CachedModel<PF,  Model<real64>>>();
    doc.render<CachedModel<PF,  Model<real32>>>();
    render_lru<3, real32>(doc);
    render_engine<MFE, 3, 0>(doc, pack<real32>(), as_pack<EnsembleType>());
}

/******************************************************************************************/

}


namespace nupack {
    void render_thermo(Document &doc) {thermo::render_mfe(doc); thermo::render_pf(doc);}
}

