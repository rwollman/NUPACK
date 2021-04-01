#include "Bind.h"
#include <nupack/execution/Local.h>
#include <nupack/types/IO.h>
#include <nupack/types/Sequence.h>
#include <nupack/types/Domain.h>
#include <nupack/types/Structure.h>
#include <nupack/common/Costs.h>

namespace nupack {

/******************************************************************************************/

void render(Document &doc, Type<Local> t) {
    doc.type(t, "core.Local");
    doc.method<0>(t, "new", rebind::construct<usize>(t));
    doc.method(t, "n_workers", &Local::n_workers);
}

/******************************************************************************************/

void render_constants(Document &doc) {
    doc.render<Base>();
    doc.render<Sequence>();
    doc.render<Strand>();
    doc.render<Domain>();
    doc.render<NamedStrand>();
    doc.render<NamedComplex>();
    doc.render<TargetStrand>();
    doc.render<Complex>();
    doc.render<PairList>();
    doc.render<Structure>();
    doc.render<TargetComplex>();

    doc.function("constants.ldexp", [](float f, std::int32_t i) {
        print(bool(i));
    });

    doc.function("constants.read_lines", [](std::string_view path) {
        std::ifstream ifs{path.data()};
        if (!ifs.good()) NUPACK_ERROR("File does not exist", path);
        vec<std::string> lines;
        for (std::string s; std::getline(ifs, s);) lines.emplace_back(s);
        std::cout << lines.front() << std::endl;
        return lines;
    });

    doc.function("core.test_matrix", [](Mat<double> const &x) {return x;});

    doc.render<SparsePairs<real>>();
    doc.function("core.sparse_pair_matrix", sparse_pair_matrix<real>);

    /// Utility functions
    doc.function("constants.dp_to_pairs", [](std::string_view s) {return io::to_pairs(s);});

    doc.function("constants.unit_evaluation_cost_table", unit_evaluation_cost_table);
    doc.function("constants.unit_evaluation_costs", unit_evaluation_costs);
    doc.function("constants.unit_subblock_cost", unit_subblock_cost);
    doc.function("constants.subblock_cost", subblock_cost<small_vec<std::size_t>>);

    doc.function("constants.trim_cxx", [](std::string s) {return trim_type_name(std::move(s), 10000);});
    doc.function("constants.rotational_symmetry", &rotational_symmetry<small_vec<uint>>);
    doc.function("constants.compute_necklaces", [](rebind::AnnotatedCallback<void, small_vec<uint>> f, uint size, uint n) {
        return compute_necklaces(small_vec<uint>(size), n, std::move(f));
    });

    // doc.function("numeric.nnls", nnls<float, std::uint16_t, float>);
    // doc.function("numeric.nnls", nnls<float, float, float>);
    // doc.function("numeric.nnls", nnls<double, double, double>);
    doc.function("constants.water_molarity", water_molarity);
    doc.function("constants.dna_salt_correction", dna_salt_correction);
    doc.object("constants.ZeroCinK", ZeroCinK);
    doc.object("constants.DefaultTemperature", DefaultTemperature);
    doc.object("constants.BoltzmannConstant", Kb);
    doc.object("constants.GitBranch", GitBranch);
    doc.object("constants.GitRevision", GitRevision);
    doc.object("constants.Version", Version);

#   define NUPACK_TMP(scope, name, val)         \
        doc.function(scope name, [] {return val;}); \
        doc.function(scope "set_" name, [](decltype(val) const &v) {val = v;})
        NUPACK_TMP("constants.", "default_parameters_path", DefaultParametersPath);
        NUPACK_TMP("constants.", "total_ram", TotalRAM);
        NUPACK_TMP("constants.", "total_cpu", TotalCPU);
#   undef NUPACK_TMP
}

void render(Document &doc, Type<AlwaysTrue> t) {doc.type(t, "constants.AlwaysTrue");}

void render(Document &doc, Type<AlwaysFalse> t) {doc.type(t, "constants.AlwaysFalse");}

void render(Document &doc, Type<True> t) {doc.type(t, "constants.TrueType");} // TODO convert to bool?

void render(Document &doc, Type<False> t) {doc.type(t, "constants.FalseType");}

}

/******************************************************************************************/

namespace rebind {

void Renderer<nupack::json>::operator()(Document &doc) const {
    using J = nupack::json;
    Type<J> t;
    doc.type(t, "core.JSON");
    doc.method(t, "new", construct<>(t));
    doc.method(t, "new", [](std::string_view s) {return J::parse(s);});
    doc.method(t, "load", [](J &j, std::string_view s) {j = J::parse(s);});
    doc.method(t, "dump", [](J const &j, unsigned indent) {return indent ? j.dump(indent) : j.dump();});

    doc.method(t, "load_binary", [](J &j, std::vector<std::uint8_t> const &s) {j = J::from_msgpack(s);});
    doc.method(t, "dump_binary", [](J const &j) {return J::to_msgpack(j);});

    doc.method(t, "load_file", [](J &j, std::string_view path) {
        std::ifstream ifs{std::string(path)};
        if (!ifs.good()) NUPACK_ERROR("invalid file", path);
        ifs >> j;
    });
}

}
