#include <nupack/design/Constraints.h>
#include <cmath>
#include <unordered_set>

namespace nupack { namespace newdesign {

using custom_csp::CompConstraint;
using custom_csp::IdentConstraint;
using custom_csp::PatternConstraint;
using custom_csp::WordConstraint;
using custom_csp::MatchConstraint;
using custom_csp::NUPACK_CS_WEAK;
using custom_csp::NUPACK_CS_STRONG;
using custom_csp::trinary;

using namespace Gecode;

using IVArray = ViewArray<Int::IntView>;
using IV = Int::IntView;

vec<int> opp_values(std::array<bool, 4> const &mask) {
    vec<int> values;
    izip(mask, [&](auto i, auto j){
        if (!j) values.emplace_back(i);
    });
    return values;
}


vec<int> nuc_values(std::array<bool, 4> const &mask) {
    vec<int> values;
    izip(mask, [&](auto i, auto j){
        if (j) values.emplace_back(i);
    });
    return values;
}


IntSet nuc_values(Base base) {
    return IntSet(nuc_values(base.mask()));
}


int random_nuc(IntVar x) {
    vec<int> choices;
    for (IntVarValues it(x); it(); ++it) {
        choices.emplace_back(it.val());
    }
    return *random_choice(choices);
}


/**
 * @brief create a string by concatenating n copies of string s for each pair (s, n) in input vector
 *
 * @param condensed vector of specifications for creating a concatenated string of
 * @return A string made by concatenating each substring in condensed its paired number of times
 */
string multiply_substrings(vec<std::pair<string, int>> const & condensed) {
    std::stringstream ss;
    for_each(condensed, [&](auto const &c) {for (auto i : range(c.second)) ss << c.first;});
    return ss.str();
}

/*****************************************************************************************************/
/*****************************************************************************************************/

/**
 * @brief A more efficient implementation of a pattern constraint propagator
 * than what can be strung together using only built-ins in Gecode
 *
 */
struct PatternProp : public Propagator {
    IVArray xs;
    // vec<vec<int>> pattern;
    Sequence pattern;

    PatternProp(Home home, IVArray _xs, Sequence _pattern)
		: Propagator(home), xs(_xs),
        // pattern(indirect_view(_pattern, [](auto x) {return nuc_values(x.mask());})) {
        pattern(_pattern) {
		xs.subscribe(home, *this, Int::PC_INT_DOM);
	}

    /* optimization for preventing posting if not a problem */
	static ExecStatus post(Home home, IVArray xs, Sequence pattern) {
        /* don't post at all if pattern can't be matched based on domains */
        for (auto i : range(xs.size())) {
            if (none_of(nuc_values(at(pattern, i).mask()), [&](auto b) {
                return xs[i].in(b);
            }))
                return ES_OK;
        }


        /* return failure if all nucleotides must match pattern */
        bool all_match = true;
        for (auto i : range(xs.size())) {
            auto poss = nuc_values(at(pattern, i).mask());

            uint count {0};
            for (auto j : poss) count += at(xs, i).in(j);
            /* there are values in the domain that aren't in the pattern at this position */
            if (count < at(xs, i).size()) {
                all_match = false;
                break;
            }
        }
        if (all_match) return ES_FAILED;

        (void) new (home) PatternProp(home, xs, pattern);
		return ES_OK;
	}

	virtual size_t dispose(Space& home) {
		xs.cancel(home, *this, Int::PC_INT_DOM);
		(void) Propagator::dispose(home);
		return sizeof(*this);
	}

	virtual void reschedule(Space& home) {
		xs.reschedule(home, *this, Int::PC_INT_DOM);
	}

	PatternProp(Space& home, PatternProp& p) : Propagator(home, p) {
		xs.update(home, p.xs);
        pattern = p.pattern;
	}

	virtual Propagator* copy(Space& home) {
		return new (home) PatternProp(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO, xs.size());
	}

    /* this is the meat of the logic */
	virtual ExecStatus propagate(Space& home, const ModEventDelta&)  {
		// if (xs.assigned())
		// 	return home.ES_SUBSUMED(*this);

        /* three lists: must match, can't match, indeterminate
        If any in can't match, then ES_SUBSUMED
        If all in must match, then ES_FAIL
        If n - 1 in must match, and 1 in indeterminate
            use .nq on all values in the pattern for the indeterminate one
        */
        vec<int> must, could;
        for (auto i : range(xs.size())) {
            auto poss = nuc_values(at(pattern, i).mask());
            auto & x = at(xs, i);
            uint count {0};
            for (auto j : poss) count += x.in(j);
            // BEEP(x, poss, count);

            if (count == 0) {
                return home.ES_SUBSUMED(*this);
            } else if (count < x.size()) {
                could.emplace_back(i);
            } else if (count == x.size()) {
                must.emplace_back(i);
            } else {
                NUPACK_BUG("impossibility in pattern constraint");
            }
        }

        // if (len(must)) BEEP(could, must);
        if (len(must) == xs.size()) {
            // BEEP(must, xs);
            return ES_FAILED;
        }

        if (len(must) == (xs.size() - 1) && len(could) == 1) {
            auto i = could[0];
            auto poss = nuc_values(at(pattern, i).mask());
            for (auto j : poss) GECODE_ME_CHECK(xs[i].nq(home, j));
        }
        return ES_FIX;
	}
};


void prevent_pattern(Home home, IntVarArgs vars, Sequence pattern) {
	NUPACK_REQUIRE(len(vars), ==, len(pattern), "mismatch in variable array and pattern size");
    GECODE_POST;
    IVArray t(home, vars);
	GECODE_ES_FAIL(PatternProp::post(home, t, pattern));
}



/*****************************************************************************************************/
/*****************************************************************************************************/


 auto remove_from_domains(IVArray window, IV allowed, WordRef ref) {
    auto const &words = ref.words();
    vec<std::array<bool, 4>> appears(len(window), Base('_').mask());

    for (auto i : range(len(words))) {
        /* only words that aren't disallowed */
        if (allowed.in(int(i))) {
            auto const &word = at(words, i);
            for (auto j : range(len(window))) for (auto val : word[j])
                appears[j][val] = true;
        }
    }

    return vmap(appears, [](auto const &app) {return opp_values(app);});
}

struct WordProp : public Propagator {
    IVArray xs;
    IV index;
    WordRef ref;

    WordProp(Home home, IVArray _xs, Int::IntView ind, WordRef _ref)
		: Propagator(home), xs(_xs), index(ind), ref(_ref) {
		xs.subscribe(home, *this, Int::PC_INT_DOM);
        index.subscribe(home, *this, Int::PC_INT_DOM);
	}


    template <class F>
	static ExecStatus post(Home home, IVArray xs, WordRef ref, F &&add_variable) {
        auto const &words = ref.words();

        /* for each word in words */
        vec<int> must, cant, could;
        izip(words, [&](auto w, auto const &word) {
            vec<int> counts;
            for (auto i : range(len(word))) {
                int count = 0;
                for (auto val : word[i]) if (xs[i].in(val)) ++count;
                if (count == 0) break;
                counts.emplace_back(count - int(xs[i].size()));
            }

            if (len(counts) < len(word)) {
                cant.emplace_back(w);
                return;
            }
            if (all_of(counts, [](auto c) {return c == 0;})) must.emplace_back(w);
            else could.emplace_back(w);
        });

        // BEEP(len(must), len(cant), len(could));

        /* rare, but if some sequence exactly matches, no need to do anything */
        if (len(must)) { return ES_OK; }
        /* also if there is only one could, this is actually a must and we can
        constrain variables further, potentially. Will cause problems if words
        have degenerate nucleotide codes */
        else if (len(could) == 1) {
            auto const &word = at(words, could[0]);
            for (auto i : range(len(xs)))
                GECODE_ME_CHECK(xs[i].eq(home, word[i][0]));
            return ES_OK;
        }
        /* no possible sequences, fail */
        else if (len(could) == 0) { return ES_FAILED; }

        /* create auxiliary variable keeping track of which words are impossible */
        IntVar temp(home, 0, len(words)-1);
        IV aux(temp);
        for (auto i : cant) GECODE_ME_CHECK(aux.nq(home, i));

        auto to_remove = remove_from_domains(xs, aux, ref);
        NUPACK_REQUIRE(len(to_remove), ==, len(xs), "length of window and potential removals must match");
        for (auto i : range(len(xs))) {
            for (auto r : at(to_remove, i)) GECODE_ME_CHECK(xs[i].nq(home, r));
        }

        add_variable(temp);
        (void) new (home) WordProp(home, xs, aux, ref);
		return ES_OK;
	}

	virtual size_t dispose(Space& home) {
		xs.cancel(home, *this, Int::PC_INT_DOM);
        index.cancel(home, *this, Int::PC_INT_DOM);
		(void) Propagator::dispose(home);
		return sizeof(*this);
	}

	virtual void reschedule(Space& home) {
		xs.reschedule(home, *this, Int::PC_INT_DOM);
        index.reschedule(home, *this, Int::PC_INT_DOM);
	}

	WordProp(Space& home, WordProp& p) : Propagator(home, p) {
		xs.update(home, p.xs);
        index.update(home, p.index);
        ref = p.ref;
	}

	virtual Propagator* copy(Space& home) {
		return new (home) WordProp(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO, xs.size() + 1);
	}

    /* this is the meat of the logic */
	virtual ExecStatus propagate(Space& home, const ModEventDelta&)  {
        // BEEP(len(*(ref.dictionary)));
        auto const &words = ref.words();

        /* if index is chosen by brancher */
        if (index.size() == 1) {
            // BEEP(index.val());
            auto const &word = at(words, index.val());
            for (auto i : range(len(xs)))
                GECODE_ME_CHECK(xs[i].eq(home, word[i][0]));
            return home.ES_SUBSUMED(*this);
        }

        /* for each word in words */
        vec<int> must, cant, could;
        izip(words, [&](auto w, auto const &word) {
            /* skip disallowed words */
            if (index.in(int(w))) {
                vec<int> counts;
                for (auto i : range(len(word))) {
                    int count = 0;
                    for (auto val : word[i]) if (xs[i].in(val)) ++count;
                    if (count == 0) break;
                    counts.emplace_back(count - int(xs[i].size()));
                }

                if (len(counts) < len(word)) {
                    cant.emplace_back(w);
                    return;
                }
                if (all_of(counts, [](auto c) {return c == 0;})) must.emplace_back(w);
                else could.emplace_back(w);
            }
        });
        // BEEP(len(must), len(cant), len(could));

        /* if one sequence must match, then this propagator is done */
        if (len(must)) { return home.ES_SUBSUMED(*this); }
        /* if only one could, it is a must */
        else if (len(could) == 1) {
            auto const &word = at(words, could[0]);
            for (auto i : range(len(xs)))
                GECODE_ME_CHECK(xs[i].eq(home, word[i][0]));
            return home.ES_SUBSUMED(*this);
        }
        /* no possible sequences, fail */
        else if (len(could) == 0) { return ES_FAILED; }
        /* remove words from index */
        for (auto i : cant) GECODE_ME_CHECK(index.nq(home, i));

        auto to_remove = remove_from_domains(xs, index, ref);
        NUPACK_REQUIRE(len(to_remove), ==, len(xs), "length of window and potential removals must match");
        for (auto i : range(len(xs))) {
            for (auto r : at(to_remove, i)) GECODE_ME_CHECK(xs[i].nq(home, r));
        }

        return ES_FIX;
	}
};

template <class F>
void constrain_window(Home home, IntVarArgs vars, WordRef word_ref, F &&f) {
	auto const &words = word_ref.words();
    for (auto &w : words) NUPACK_REQUIRE(len(vars), ==, len(w), w, "A word is a different length than the window");
    GECODE_POST;
    IVArray t(home, vars);
	GECODE_ES_FAIL(WordProp::post(home, t, word_ref, fw<F>(f)));
}

/*****************************************************************************************************/
/*****************************************************************************************************/

NucSpace::NucSpace(Sequence const &domains) {
    IntVarArgs temp;
    for (auto d : domains) {
        temp << IntVar(*this, nuc_values(d));
    }
    nucs = IntVarArray(*this, temp);
}


NucSpace::NucSpace(NucSpace &other) : Space(other) {
    nucs.update(*this, other.nucs);
    ref = other.ref;
    // ref.update(*this, other.ref);
    extras.update(*this, other.extras);
}


NucSpace * NucSpace::copy() {
    return new NucSpace(*this);
}


NucSpace * NucSpace::cast_clone() {
    status();
    return static_cast<NucSpace *>(this->clone());
}



void NucSpace::match_constraint(int i, int j) {
    rel(*this, nucs[i] == nucs[j]);
}


void NucSpace::complementarity_constraint(int i, int j, bool wobble) {
    if (!wobble) {
        rel(*this, nucs[i] == 3 - nucs[j]);
    } else {
        BoolVarArray temp(*this, 4, 0, 1);
        rel(*this, nucs[i], IRT_EQ, 0, pmi(temp[0]));
        rel(*this, nucs[j], IRT_EQ, 3, imp(temp[0]));

        rel(*this, nucs[i], IRT_EQ, 1, pmi(temp[1]));
        rel(*this, nucs[j], IRT_EQ, 2, imp(temp[1]));

        rel(*this, nucs[i], IRT_EQ, 2, pmi(temp[2]));
        dom(*this, nucs[j], nuc_values(Base('Y')), imp(temp[2]));

        rel(*this, nucs[i], IRT_EQ, 3, pmi(temp[3]));
        dom(*this, nucs[j], nuc_values(Base('R')), imp(temp[3]));


        // rel(*this, BOT_XOR, temp, 1);
        // rel(*this, nucs[i], )

        // BoolVarArgs temp;
        // temp << expr(*this, (nucs[i] == 0) && (nucs[j] == 3)); /* A */
        // temp << expr(*this, (nucs[i] == 1) && (nucs[j] == 2)); /* C */
        // temp << expr(*this, (nucs[i] == 2) && (singleton(nucs[j]) <= IntSet({1, 3}))); /* G */
        // temp << expr(*this, (nucs[i] == 3) && (singleton(nucs[j]) <= IntSet({0, 2}))); /* T */
        // rel(*this, BOT_OR, temp, 1);
    }
}


void NucSpace::pattern_constraint(vec<int> const &window, Sequence const &pattern) {
    uint n = len(pattern);
    uint end = len(window) - n + 1;
    for (auto i : range(end)) {
        auto inds = view(window, i, i+n);
        IntVarArgs temp;
        for (auto j : inds) temp << nucs[j];

        prevent_pattern(*this, temp, pattern);
        // BoolVarArgs temp;

        // /* for each nucleotide in subwindow, check if its values is in the domain of the pattern nucleotide */
        // izip(inds, pattern, [&](auto ind, auto i, auto d) {
        //     BoolVar b(*this, 0, 1);
        //     dom(*this, nucs[i], nuc_values(d), b);
        //     temp << b;
        //     // temp << expr(*this, singleton(nucs[i]) <= nuc_values(d));
        // });

        // /* add the constraint */
        // rel(*this, BOT_AND, temp, 0);
    }
}


void NucSpace::diversity_constraint(vec<int> const &window, int wordsize, int mintypes) {
    uint end = len(window) - wordsize + 1;
    for (auto i : range(end)) {
        auto inds = view(window, i, i+wordsize);
        IntVarArgs temp;
        for (auto j : inds) temp << nucs[j];

        nvalues(*this, temp, IRT_GQ, mintypes);
    }
}


// void NucSpace::word_constraint(vec<int> const &window, vec<Sequence> const &words) {
void NucSpace::word_constraint(vec<int> const &window, WordRef word_ref) {
    IntVarArgs temp;
    for (auto i : window) temp << nucs[i];
    constrain_window(*this, temp, word_ref, [&](auto var) {
        IntVarArgs temp(this->extras);
        temp << var;
        this->extras = IntVarArray(*this, temp);
        // BEEP(len(this->extras));
    });

    // int n = len(window);
    // BoolVarArgs word_matches;
    // for (auto const &w : word_ref.words()) {
    //     BoolVarArgs temp;
    //     zip(window, w, [&](auto i, auto d) {
    //         temp << expr(*this, singleton(nucs[i]) == nuc_values(d));
    //     });
    //     word_matches << expr(*this, sum(temp) == n);
    // }
    // rel(*this, BOT_OR, word_matches, 1);
}


void NucSpace::similarity_constraint(vec<int> const &window, WordRef word_ref, std::pair<int, int> limits) {
    // IntVarArgs var_window;
    // for (auto i : window) var_window << nucs[i];

    // similar(*this, var_window, word_ref, limits);

    BoolVarArgs temp;
    auto const &reference = word_ref.word();
    zip(window, reference, [&](auto i, auto d) {
        // BoolVar b(*this, 0, 1);
        BoolVarArgs indiv_match;
        for (auto v : d) {
            auto match = expr(*this, nucs[i] == v);
            indiv_match << match;
        }
        // rel(*this, BOT_OR, indiv_match, 1, b);
        temp << expr(*this, sum(indiv_match) >= 1);
        // dom(*this, nucs[i], nuc_values(d), b);
        // temp << expr(*this, singleton(nucs[i]) <= nuc_values(d));
    });

    auto [lower, upper] = limits;
    if (lower == upper) rel(*this, sum(temp) == upper);
    else {
        rel(*this, sum(temp) <= upper);
        rel(*this, sum(temp) >= lower);
    }
}



void NucSpace::force(int i, int val) {
    rel(*this, nucs[i] == val);
}


void NucSpace::disallow(int i, int val) {
    rel(*this, nucs[i] != val);
}


void NucSpace::add_reference(Sequence const &reference) {
    // IntVarArgs temp;
    for (auto d : reference) {
        // int value {d.value};
        NUPACK_REQUIRE(d.value, <, 4, "Degenerate reference sequence");
        // temp << IntVar(*this, IntSet(vec<int>{value}));
    }
    ref = &reference;
    // ref = IntVarArray(*this, temp);
}


auto select_close = [](auto const &home, auto x, auto i) {
    NucSpace const &space = static_cast<NucSpace const &>(home);
    auto const &ref = *(space.ref);
    int ret {-1};
    if (i < ref.size()) {
        auto r = ref[i].value;
        if (x.in(r)) ret = r;
    }
    if (ret == -1) ret = random_nuc(x);
    return ret;
};

auto val_select() {
    return INT_VAL(select_close);
};



void NucSpace::default_brancher() {
    Rnd r(random_range(1, 1e10));
    branch(*this, extras, INT_VAR_NONE(), INT_VAL_RND(r));

    IntAFC afc(*this, nucs, 0.99);
    // IntAction act(*this, nucs, 0.99);
    // branch(*this, nucs, tiebreak(INT_VAR_ACTION_MAX(act), INT_VAR_SIZE_MIN(), INT_VAR_RND(r)), INT_VAL_RND(r));
    branch(*this, nucs, tiebreak(INT_VAR_AFC_MAX(afc), INT_VAR_RND(r)), INT_VAL_RND(r));
    // branch(*this, nucs, INT_VAR_RND(r), INT_VAL_RND(r));
    // branch_rnd_afc(*this, nucs);
    // branch(*this, nucs, tiebreak(INT_VAR_SIZE_MIN(), INT_VAR_RND(r)), INT_VAL_RND(r));
}




void NucSpace::reference_brancher() {
    Rnd r(random_range(1, 1e10));
    IntAFC afc(*this, nucs, 0.99);
    // IntAction act(*this, nucs, 0.99);
    // branch(*this, nucs, tiebreak(INT_VAR_ACTION_MAX(act), INT_VAR_SIZE_MIN(), INT_VAR_RND(r)), val_select());
    // branch(*this, nucs, INT_VAR_NONE(), val_select());
    // branch(*this, nucs, INT_VAR_RND(r), val_select(ref));
    // branch(*this, nucs, tiebreak(INT_VAR_AFC_MAX(afc), INT_VAR_SIZE_MIN(), INT_VAR_RND(r)), val_select());
    branch(*this, nucs, tiebreak(INT_VAR_SIZE_MIN(), INT_VAR_RND(r)), val_select());
    // branch_rnd_afc(*this, nucs);
    // heuristic_branch(*this, nucs, &state);
}

void NucSpace::cheap_reference_brancher() {
    branch(*this, nucs, INT_VAR_NONE(), val_select());
}


NucSpace::operator Sequence() const {
    return vmap<Sequence>(nucs, [](auto i) {return Base::from_index(i.val());});
}



Constraints::Constraints(Sequence const &domains) : initial(std::make_unique<NucSpace>(domains)), handler() {
    for_each(domains, [&](auto const &n) {handler.add_nucleotide_variable(n);});
};





void Constraints::match_constraint(int i, int j) {
    /* new */
    initial->match_constraint(i, j);

    /* old */
    handler.add_constraint(IdentConstraint(i, j));
}


void Constraints::complementarity_constraint(int i, int j, bool wobble) {
    /* new */
    initial->complementarity_constraint(i, j, wobble);

    /* old */
    if (wobble) {
        handler.add_constraint(CompConstraint(i, j, NUPACK_CS_WEAK));
    } else {
        handler.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
    }
}


void Constraints::pattern_constraint(vec<int> const &window, Sequence const &pattern) {
    /* do nothing if pattern won't be a problem */
    if (len(pattern) > len(window)) return;

    /* new */
    initial->pattern_constraint(window, pattern);

    /* old */
    auto poss = handler.get_possible_nucleotides();
    handler.add_constraint(PatternConstraint(window, string(pattern), poss));
}


void Constraints::diversity_constraint(vec<int> const &window, int wordsize, int mintypes) {
    /* new */
    initial->diversity_constraint(window, wordsize, mintypes);

    /* old */
    static const std::map<int, vec<string>> diversity_levels {
        {2, {"A", "C", "T", "G"}},
        {3, {"W", "S", "M", "K", "Y", "R"}},
        {4, {"V", "H", "D", "B"}}
    };

    if (mintypes > 1 && mintypes <= 4) {
        auto poss = handler.get_possible_nucleotides();

        for (auto s : diversity_levels.at(mintypes)) {
            auto pattern = multiply_substrings({{s, wordsize}});
            handler.add_constraint(PatternConstraint(window, pattern, poss));
        }
    }
}


void Constraints::word_constraint(vec<int> const &window, vec<Sequence> const &words) {
    // for (auto const &w: words) NUPACK_REQUIRE(len(window), ==, len(w), "word/window size mismatch");
    // initial->word_constraint(window, words);
    /* new */
    DictWords int_words;
    for (auto const &word: words) {
        int_words.emplace_back(vmap(word, [](auto w) {return nuc_values(w.mask());}));
    }

    dictionary->emplace_back(int_words);
    initial->word_constraint(window, WordRef{dictionary.get(), uint(dictionary->size() - 1)});

    /* old */
    auto supp_var = handler.add_variable(vec<trinary>(len(words), true));
    auto old_words = vmap(words, [](auto const &x) {return string(x);});
    handler.add_constraint(WordConstraint(window, old_words, supp_var));
    ++num_extra_vars;
}


void Constraints::similarity_constraint(vec<int> const &window, Sequence const &reference, std::pair<real, real> limits) {
    NUPACK_REQUIRE(len(window), ==, len(reference), "window/reference size mismatch");
    NUPACK_REQUIRE(limits.first, >=, 0, "lower limit invalid");
    NUPACK_REQUIRE(limits.first, <, limits.second, "limits out of order");
    NUPACK_REQUIRE(limits.first, <=, 1, "upper limit invalid");

    uint n = len(window);
    int lower = std::ceil(limits.first * n);
    int upper = std::floor(limits.second * n);
    NUPACK_REQUIRE(lower, <=, upper, "discretized limits are incompatible");
    DictWords one;
    one.emplace_back(vmap(reference, [](auto w) {return nuc_values(w.mask());}));
    dictionary->emplace_back(one);

    /* new */
    initial->similarity_constraint(window, WordRef{dictionary.get(), uint(dictionary->size() - 1)}, {lower, upper});

    /* old */
    handler.add_constraint(MatchConstraint(window, string(reference), {limits.first}, {limits.second}));
}






Optional<Sequence> Constraints::initial_sequence() {
    /* try new */
    auto new_res = new_initial_sequence();
    if (std::holds_alternative<Sequence>(new_res))
        return std::get<Sequence>(new_res);
    else {
        /* no possible sequences */
        if (!std::get<bool>(new_res)) return optional();
    }

    /* try old */
    auto old_res = old_initial_sequence();
    if (old_res) return old_res.value();

    return optional();
}

Optional<Sequence> Constraints::make_mutation(Sequence const &seq, vec<int> positions) {
    auto cur_seq = seq;

    for (auto const &pos : positions) {
        // The way this is currently done is to try the new solver and stop after some time limit.
        // This makes it non-deterministic.
        // It would be better to adopt a non-time based metric or just make a better constraint solver
        // that doesn't have to use this switching behavior.
        // It does seem that the new constraint solver is really slow (maybe even never finishes)
        // On the other hand it's faster in most cases I think.
        // So the approach now is just to set msec_cutoff to 0 for deterministic designs.

        /* try new */
        if (msec_cutoff > 0) { // disable time based switching if deterministic seed.
            auto new_res = make_new_mutation(cur_seq, pos);
            if (std::holds_alternative<Sequence>(new_res)) {
                cur_seq = std::get<Sequence>(new_res);
                continue;
            } else {
                /* no possible sequences */
                if (!std::get<bool>(new_res)) continue;
            }
        }

        /* try old */
        auto old_res = make_old_mutation(cur_seq, pos);
        if (old_res) cur_seq = old_res.value();
    }

    /* at least one mutation worked */
    if (cur_seq != seq) return cur_seq;

    return optional();
}

Variant<Sequence, bool> Constraints::make_new_mutation(Sequence const &cur_seq, int pos) {
    auto opts = search_options();

    std::unique_ptr<NucSpace> space_ptr(initial->cast_clone());
    NucSpace &space = *space_ptr;

    /* add reference and disallow current variable-value assignment */
    space.add_reference(cur_seq);

    /* change space so that nucleotide must be mutated */
    space.disallow(pos, at(cur_seq, pos).value);

    /* set branching based on whether time-based stop is present */
    space.cheap_reference_brancher();

    DFS<NucSpace> searcher(&space, opts);
    std::unique_ptr<NucSpace> result(searcher.next());

    if (result) { return Sequence(*result); }
    else if (searcher.stopped()) { return true; }
    else { return false; }
}

Optional<Sequence> Constraints::make_old_mutation(Sequence const &cur_seq, int pos) {
    Optional<Sequence> ret;

    auto time = time_it([&] {
        vec<int> cur_vars{view(cur_seq)};
        /* correct length when window constraints exist */
        for (auto i : range(num_extra_vars)) cur_vars.emplace_back(0);

        auto vars = handler.make_mutation({pos}, cur_vars);

        if (vars[0] == -1) {
            ret = optional();
            return;
        }

        if (len(vars) > sequence_length()) vars.resize(sequence_length());
        ret = vmap<Sequence>(vars, [](auto i) {return Base::from_index(i);});
    });

    update_cutoff(time);
    return ret;
}

Variant<Sequence, bool> Constraints::new_initial_sequence() {
    std::unique_ptr<NucSpace> space_ptr(initial->cast_clone());
    NucSpace &space = *space_ptr;

    /* add empty reference to avoid segfault */
    Sequence seq;
    space.add_reference(seq);
    space.default_brancher();

    auto opts = search_options();

    DFS<NucSpace> searcher(&space, opts);
    std::unique_ptr<NucSpace> result(searcher.next());

    if (result) { return Sequence(*result); }
    else if (searcher.stopped()) { return true; }
    else { return false; }
}

Optional<Sequence> Constraints::old_initial_sequence() {
    Optional<Sequence> ret;

    auto time = time_it([&]{
        auto vars = handler.init_random();
        if (len(vars) > sequence_length()) vars.resize(sequence_length());
        ret = vmap<Sequence>(vars, [](auto i) {return Base::from_index(i);});
    });

    update_cutoff(time);
    return ret;
}

int Constraints::sequence_length() const {
    return len(initial->nucs);
}


Gecode::Search::Options Constraints::search_options() const {
    Search::Options options;
    // options.c_d = len(initial->nucs) / 10;
    options.c_d = len(initial->nucs) * 2;
    // options.a_d = 10;
    options.stop = Search::Stop::time(msec_cutoff);
    return options;
}


void Constraints::update_cutoff(real time) {
    old_mut_time.add_value(time);
    auto avg = old_mut_time.average;
    if (msec_cutoff > 0) msec_cutoff = int(avg * 100);
}

real RunningAverage::add_value(real val) {
    average = (val + (average * count)) / (count + 1);
    ++count;
    return average;
}

}}
