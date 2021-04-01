#include <nupack/thermo/CachedModel.h>
#include <nupack/thermo/Subopt.h>
#include <nupack/model/Model.h>
#include <nupack/state/State.h>
#include <nupack/types/Structure.h>

namespace nupack::thermo {

/******************************************************************************************/

/******************************************************************************************/

// Return unique secondary structures with their structure energy and lowest stack energy
// Free energies are given wrt indistinguishable strands in the current API.
std::map<Structure, std::pair<real, real>> unique_subopt(vec<std::pair<PairList, real>> v, Complex const &c, Model<float> const &model) {
    std::map<Structure, std::pair<real, real>> out;
    auto const nicks = c.nicks();

    if (model.ensemble == Ensemble::stacking) {
        auto const sys = std::make_shared<System const>(c.strands());
        for (auto &p : v) {
            auto [it, put] = out.try_emplace(Structure(std::move(p.first), nicks));
            it->second.first = structure_energy(sys, it->first, model);
            min_eq(it->second.second, p.second);
        }
    } else {
        for (auto &p : v) out.try_emplace(Structure(std::move(p.first), nicks), p.second, p.second);
    }

    // Do NOT include rotational symmetry now

    // if (auto const sym = rotational_symmetry(c.strands()); sym != 1) {
    //     real dE = std::log(real(sym)) / model.beta;
    //     for (auto &p : out) {
    //         p.second.first += dE;
    //         p.second.second += dE;
    //     }
    // }

    return out;
}

/******************************************************************************************/

}
