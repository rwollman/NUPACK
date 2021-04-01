#pragma once
#include <cmath>

namespace nupack {

/******************************************************************************************/

template <class V = vec<int>>
vec<vec<int>> get_masks(V const & lengths) {
    vec<int> stackable;
    izip(lengths, [&stackable](auto i, auto const & l) {if (l == 2) stackable.emplace_back(i);});

    int num = std::pow(2, len(stackable));
    vec<vec<int>> masks(num);
    // https://stackoverflow.com/questions/2686542/converting-integer-to-a-bit-representation
    // find all unconstrained stackings
    for (auto x : range(num)) {
        auto & ret = masks[x];
        while (x) {
            ret.push_back(int(x & 1));
            x >>= 1;
        }
        ret.resize(len(stackable), 0);
    }

    // spread out maps onto Full length vector.
    auto full_masks = vmap(masks, [&](auto const & mask) {
        vec<int> full_mask(len(lengths));
        zip(mask, stackable, [&](auto m, auto s) {if (m) full_mask[s] = 1;});
        return full_mask;
    });

    return full_masks;
}

/******************************************************************************************/

// remove stack configurations with consecutive stacks => conflicting stacks on a pair.
template <class V = vec<vec<int>>>
void filter_masks(V & masks) {
    erase_if(masks, [](auto const & c) {
        bool check = false;
        zip(view(c, 0, len(c) - 1), view(c, 1, len(c)), [&](auto l, auto r) {
            if (l && r) check = true;
        });
        return check;
    });
}

/******************************************************************************************/
/*              Full enumerated partition function with optional counting                 */
/******************************************************************************************/

// compute the number of dangle configurations for each coaxial stacking state
/// partition function over the set of masks
template <class Rig, class V = vec<vec<int>>, class S, class M>
real mask_pf(V const & masks, S seqs, M const & m, bool circ, real left, real right) {

    auto lengths = vmap(seqs, len);
    if (!circ) lengths = vmap(view(seqs, 1, len(seqs) - 1), len);

    if (circ) {
        seqs.insert(begin_of(seqs), copy(back(seqs)));
        seqs.push_back(seqs[1]);
    }

    real total = Rig::zero();
    real const one = Rig::one();

    // sum over each mask
    for (auto mask : masks) {

        // multiply in outermost dangles if outermost stacks are inactive
        // logic reduces to outside = 1 for multiloops
        real product = Rig::times(mask.front() ? one : left, mask.back() ? one : right);

        // default for exterior loop
        int first = 0, last = 0;
        // circular; multiloop
        if (circ) {first = mask.front(); last = mask.back();}
        // for creating shifted views
        mask.insert(begin_of(mask), last);
        mask.insert(end_of(mask), first);

        zip(range(len(lengths)).shift(1), lengths, view(mask, 0, len(mask) - 2), view(mask, 2, len(mask)),
            view(mask, 1, len(mask) - 1), [&](auto i, auto length, auto l, auto r, auto c) {
            if (c) {
                product = Rig::times(product,
                    Rig::boltz(m.beta, m.coaxial_stack_energy(back(seqs[i - 1]), front(seqs[i]),
                                                              back(seqs[i]), front(seqs[i + 1]))));
            }
            if (length >= 3) {
                auto left = [&]{
                    return Rig::boltz(m.beta, m.dG(dangle5, back(seqs[i-1]), front(seqs[i]), front(seqs[i], 1)));
                };
                auto right = [&]{
                    return Rig::boltz(m.beta, m.dG(dangle3, back_index(seqs[i], 1), back(seqs[i]), front(seqs[i+1])));
                };

                if (!l && r) {
                    product = Rig::times(product, Rig::plus(one, left()));
                } else if (l && !r) {
                    product = Rig::times(product, Rig::plus(one, right()));
                } else if (l + r == 0) {
                    if (length == 3) {
                        product = Rig::times(product, Rig::plus(one, left(), right()));
                    }
                    else if (length > 3) {
                        product = Rig::times(product, Rig::plus(one, left()), Rig::plus(one, right()));
                    }
                }
            }
        });
        total = Rig::plus(total, product);
    };
    return total;
}

/******************************************************************************************/

/// stacking partition function of exterior loop
template <class Rig, class V, class M>
real exterior_loop_stack_sum(V seqs, M const & m) {
    if (len(seqs) == 1) return Rig::one();

    // pre-underscore first; post-underscore last.
    std::rotate(begin_of(seqs), max_element(seqs), end_of(seqs));

    real left = Rig::one(), right = Rig::one();

    // dangle states for 5' end of loop
    auto first = seqs.front();
    if (len(first) > 2) {
        auto & second = seqs[1];
        left = Rig::plus(Rig::one(), Rig::boltz(m.beta, m.dG(dangle3, back_index(first, 1), back(first), second[0])));
    }
    // dangle states for 3' end of loop
    auto last = seqs.back();
    if (len(last) > 2) {
        auto & penultimate = back_index(seqs, 1);
        right = Rig::plus(Rig::one(), Rig::boltz(m.beta, m.dG(dangle5, back(penultimate), last[0], last[1])));
    }

    if (len(seqs) == 2) { // only dangles
        return Rig::times(left, right);
    } else if (len(seqs) > 2) {
        auto mask_seqs = view(seqs, 1, len(seqs) - 1);
        auto lengths = vmap(mask_seqs, len);
        auto masks = get_masks(lengths);
        filter_masks(masks);
        return mask_pf<Rig>(masks, seqs, m, false, left, right);
    }
    return Rig::one();
}

/******************************************************************************************/

/// stacking partition function of multi loop
template <class Rig, class V, class M>
real multiloop_stack_sum(V const &seqs, M const & m) {
    auto lengths = vmap(seqs, len);
    auto masks = get_masks(lengths);

    // so that first and last stacks are consecutive for filtering
    for (auto & mask : masks) mask.push_back(mask[0]);
    filter_masks(masks);
    // remove auxiliary repeated first base
    for (auto & mask : masks) mask.pop_back();

    return mask_pf<Rig>(masks, seqs, m, true, Rig::one(), Rig::one());
}

/******************************************************************************************/

/// partition function of any loop
template <class Rig, class O, class M>
real loop_stack_sum(O const & o, M const & m) {
    if (!o.exterior() && len(o.seqs) < 3) { // hairpin, interior, bulge, stack
        return Rig::boltz(m.beta, m.loop_energy(o.sequences(), o.nick()));
    } else if (o.exterior()) {
        return exterior_loop_stack_sum<Rig>(o.sequences(), m);
    } else {
        return Rig::times(Rig::boltz(m.beta, m.linear_multi_energy(o.sequences())),
                          multiloop_stack_sum<Rig>(o.sequences(), m));
    }
}

/******************************************************************************************/

}
