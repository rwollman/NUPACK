#include "Math.h"

namespace nupack {

namespace concentration {

/**************************************************************************************/

void render(Document &doc, Type<Options> t);

/**************************************************************************************/

template <class V>
void render(Document &doc, Type<Output<V>> t) {
    doc.type(t, "concentration.Output");
    render_public(doc, t);
}

/**************************************************************************************/

rebind::Variable response(std::type_index, Options const &v) {return to_dictionary(v);}

/**************************************************************************************/

std::optional<Options> request(Type<Options>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto z = r.request<rebind::Dictionary>()) return from_dictionary<Options>(std::move(*z), msg);
    return msg.error("Not a dictionary-like type");
}

/**************************************************************************************/

void render(Document &doc, Type<Options> t) {
    doc.type(t, "concentration.Options");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "set_method", [](Options &o, uint n) {NUPACK_REQUIRE(n, <, 6); o.method = static_cast<Method>(n);});
    render_public(doc, t);
}

}

/**************************************************************************************/

void render_math(Document &doc) {
    doc.function("concentration.solve", &concentration::equilibrate);
    doc.function("concentration.solve_complexes", &concentration::solve_complexes);
}

/**************************************************************************************/

}
