#include <nupack/model/Model.h>
#include <nupack/model/ParameterSet.h>
#include <nupack/common/Error.h>
#include <nupack/common/Runtime.h>
#include <nupack/state/State.h>
#include <nupack/reflect/Serialize.h>
#include <fstream>

namespace nupack {

/******************************************************************************************/

ParameterFile::ParameterFile(string p) : path(p) {
    if (path == "RNA" || path == "rna") path = "rna06";
    if (path == "DNA" || path == "dna") path = "dna04";
    if (path.size() < 5 || path.substr(path.size() - 5) != ".json") path += ".json";
}

/******************************************************************************************/

json ParameterFile::open() const {
    string name = path;
    // BEEP(DefaultParametersPath);
    if (!path_exists(name)) {
        name = path_join(DefaultParametersPath, path);
    }
    if (!path_exists(name)) {
        auto s = get_env("NUPACKHOME");
        if (!s.empty()) name = path_join(path_join(s, "parameters"), path);
    }

    std::ifstream file(name);
    if (!file.good()) {
        vec<string> directories = {".", DefaultParametersPath, "$NUPACKHOME/parameters"};
        NUPACK_ERROR("failed to open parameter file ", path, directories);
    }
    json j;
    file >> j;
    return j;
}

/******************************************************************************************/

std::array<char const *, 5> EnsembleNames = {"nostacking", "stacking",  "min", "all", "none"};

Ensemble as_ensemble(string_view s) {
    if (auto it = find(EnsembleNames, s); it != EnsembleNames.end())
        return Ensemble(it - EnsembleNames.begin());
    NUPACK_ERROR("invalid ensemble type", s);
}

/******************************************************************************************/

int find_loop_structure_nick(Complex const &c, PairList const &p) {
    NUPACK_REQUIRE(len(c), ==, len(p));
    int nick = -1;
    auto const n_pairs = p.n_pairs();

    for (auto const &s : c.views()) NUPACK_REQUIRE(len(s), >, 0, "Strands should have nonzero length");

    // find a nick ...
    if (c.size() == 2 && c.n_strands() == 2) {
        nick = 0; // horrific edge case.
    } else {
        izip(c.positions, [&](auto i, auto n) { // n is first index of the strand after i
            if (p[n-1] != n % len(c))
                nick = (i+1) % c.n_strands();
        });
        if (n_pairs + int(nick != -1) != c.n_strands()) {
            NUPACK_ERROR("Incorrect number of base pairs for a loop secondary structure", n_pairs, c.n_strands(), nick);
        }
    }


    return nick;
}

/******************************************************************************************/

// Insert null bases if there is a nick and they don't exist
SequenceList complex_to_loop(Complex const &c, int nick) {
    auto v = vmap<SequenceList>(c.views(), [](auto const &s) {return Sequence(s);});
    NUPACK_REQUIRE(nick, >=, -2);
    NUPACK_REQUIRE(nick, <, int(len(v)));
    for (auto const &s : v)
        NUPACK_ASSERT(!s.empty(), "Loop contains an empty sequence");

    if (nick >= 0) {
        auto &s5 = v[nick];
        if (front(s5) != Base('_')) s5.insert(s5.begin(), Base('_'));
        auto &s3 = nick ? v[nick - 1] : v.back();
        if (back(s3) != Base('_')) s3.push_back(Base('_'));
    }

    for (auto const &s : v)
        NUPACK_REQUIRE(s.size(), >=, 2, "Loop sequence does not contain enough nucleotides", s, v);

    return v;
}

/******************************************************************************************/

}
