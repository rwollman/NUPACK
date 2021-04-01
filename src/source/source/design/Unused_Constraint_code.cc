/*****************************************************************************************************/
/*****************************************************************************************************/

/* vec implementation */
struct Mask {
    vec<bool> mask;
    int num_true;

    Mask() = default;
    Mask(int length) : mask(length, false), num_true{0} {}

    void emplace(int mark) {
        if (!mask[mark]) {
            mask[mark] = true;
            ++num_true;
        }
    }

    bool contains(int i) const { return mask[i]; }

    auto size() const { return static_cast<size_t>(num_true); };
};

// using Mask = std::set<int>;

struct SimilarityProp : public Propagator {
    IVArray xs;
    Mask cant;
    Mask must;
    WordRef ref;
    std::pair<int, int> lims;

    SimilarityProp(Home home, IVArray _xs, Mask _cant, Mask _must,
            WordRef _ref, std::pair<int, int> _lims)
		: Propagator(home), xs(_xs), cant(_cant), must(_must), ref(_ref), lims(_lims) {
		xs.subscribe(home, *this, Int::PC_INT_DOM);
	}


	static ExecStatus post(Home home, IVArray xs, WordRef ref, std::pair<int, int> lims) {
        auto const & word = ref.word();
        auto n = len(word);
        Mask cant(n), must(n);

        for (auto i : range(len(xs))) {
            int count = 0;
            for (auto val : word[i]) if (xs[i].in(val)) ++count;

            if (count == 0) cant.emplace(i);
            else if (count == len(xs[i])) must.emplace(i);
        }

        /* too many or too few matching */
        if (len(must) > lims.second || (n - len(cant)) < lims.first) return ES_FAILED;

        /* exactly at upper limit, no more allowed */
        if (len(must) == lims.second) {
            for (auto i : range(len(xs))) {
                bool undetermined = !must.contains(i) && !cant.contains(i);
                if (undetermined) {
                    for (auto val : word[i]) GECODE_ME_CHECK(xs[i].nq(home, int(val)));
                }
            }
            return ES_OK;
        }

        /* exactly at upper limit of failure, all remaining must match */
        if (n - len(cant) == lims.first) {
            std::array<bool, 4> mask;
            for (auto i : range(len(xs))) {
                bool undetermined = !must.contains(i) && !cant.contains(i);
                if (undetermined) {
                    /* reset */
                    for (auto &v : mask) v = true;
                    /* remove reference from mask */
                    for (auto val : word[i]) mask[val] = false;
                    /* disallow non-matching */
                    for (auto j : range(len(mask)))
                        if (mask[j]) GECODE_ME_CHECK(xs[i].nq(home, int(j)));
                }
            }
            return ES_OK;
        }

        (void) new (home) SimilarityProp(home, xs, cant, must, ref, lims);
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

	SimilarityProp(Space& home, SimilarityProp& p) : Propagator(home, p) {
		xs.update(home, p.xs);
        cant = p.cant;
        must = p.must;
        ref = p.ref;
        lims = p.lims;
	}

	virtual Propagator* copy(Space& home) {
		return new (home) SimilarityProp(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO, xs.size() + 2);
	}

    bool undetermined(int i) const {
        return !must.contains(i) && !cant.contains(i);
    }

	virtual ExecStatus propagate(Space& home, const ModEventDelta&)  {
        auto const & word = ref.word();
        auto n = len(word);

        for (auto i : range(len(xs))) {
            if (undetermined(i)) {
                int count = 0;
                for (auto val : word[i]) if (xs[i].in(val)) ++count;

                if (count == 0) cant.emplace(i);
                else if (count == len(xs[i])) must.emplace(i);
            }
        }

        /* too many or too few matching */
        if (len(must) > lims.second || (n - len(cant)) < lims.first) return ES_FAILED;

        /* exactly at upper limit, no more allowed */
        if (len(must) == lims.second) {
            for (auto i : range(len(xs))) {
                if (undetermined(i)) {
                    for (auto val : word[i]) GECODE_ME_CHECK(xs[i].nq(home, int(val)));
                }
            }
            return home.ES_SUBSUMED(*this);
        }

        /* exactly at upper limit of failure, all remaining must match */
        if (n - len(cant) == lims.first) {
            std::array<bool, 4> mask;
            for (auto i : range(len(xs))) {
                if (undetermined(i)) {
                    /* reset */
                    for (auto &v : mask) v = true;
                    /* remove reference from mask */
                    for (auto val : word[i]) mask[val] = false;
                    /* disallow non-matching */
                    for (auto j : range(len(mask)))
                        if (mask[j]) GECODE_ME_CHECK(xs[i].nq(home, int(j)));
                }
            }
            return home.ES_SUBSUMED(*this);
        }
        return ES_FIX;
	}
};

void similar(Home home, IntVarArgs vars, WordRef word_ref, std::pair<int, int> lims) {
    GECODE_POST;
    IVArray t(home, vars);
	GECODE_ES_FAIL(SimilarityProp::post(home, t, word_ref, lims));
}

/*****************************************************************************************************/
/*****************************************************************************************************/

/*****************************************************************************************************/
/*****************************************************************************************************/

/* Get this to implement a weighted random branching with AFC as weight */

class RandomAFC : public Brancher {
protected:
    IVArray x;
    // choice definition
    class PosVal : public Choice {
    public:
        int pos; int val;
        PosVal(const RandomAFC& b, int p, int v)
            : Choice(b,2), pos(p), val(v) {}
        virtual void archive(Archive& e) const {
            Choice::archive(e);
            e << pos << val;
        }
    };
public:
    RandomAFC(Home home, IVArray& x0)
        : Brancher(home), x(x0) {}
    static void post(Home home, IVArray& x) {
        (void) new (home) RandomAFC(home,x);
    }
    virtual size_t dispose(Space& home) {
        (void) Brancher::dispose(home);
        return sizeof(*this);
    }
    RandomAFC(Space& home, RandomAFC& b)
        : Brancher(home,b) {
        x.update(home,b.x);
    }
    virtual Brancher* copy(Space& home) {
        return new (home) RandomAFC(home,*this);
    }
    // status
    virtual bool status(const Space& home) const {
        for (int i=0; i<x.size(); i++)
            if (!x[i].assigned())
                return true;
        return false;
    }
    // choice
    virtual Choice* choice(Space& home) {
        vec<real> pos_vars;
        pos_vars.reserve(len(x));
        for (auto i: range(len(x))) {
            if (x[i].assigned()) pos_vars.emplace_back(0);
            else {
                // pos_vars.emplace_back(5-x[i].size());
                pos_vars.emplace_back(x[i].size());
                // pos_vars.emplace_back(x[i].afc());
            }
        }
        auto it = random_choice(pos_vars);
        auto i = it - begin_of(pos_vars);
        auto val = select_close(home, x[i], i);
        return new PosVal(*this, i, val);

        GECODE_NEVER;
        return NULL;
    }
    virtual Choice* choice(const Space&, Archive& e) {
        int pos, val;
        e >> pos >> val;
        return new PosVal(*this, pos, val);
    }
    // commit
    virtual ExecStatus commit(Space& home,
                                                        const Choice& c,
                                                        unsigned int a) {
        const PosVal& pv = static_cast<const PosVal&>(c);
        int pos=pv.pos, val=pv.val;
        if (a == 0)
            return me_failed(x[pos].eq(home,val)) ? ES_FAILED : ES_OK;
        else
            return me_failed(x[pos].nq(home,val)) ? ES_FAILED : ES_OK;
    }
    // print
    virtual void print(const Space& home, const Choice& c,
                                         unsigned int a,
                                         std::ostream& o) const {
        const PosVal& pv = static_cast<const PosVal&>(c);
        int pos=pv.pos, val=pv.val;
        if (a == 0)
            o << "x[" << pos << "] = " << val;
        else
            o << "x[" << pos << "] != " << val;
    }
};

void branch_rnd_afc(Home home, const IntVarArgs& x) {
    if (home.failed()) return;
    IVArray y(home,x);
    RandomAFC::post(home,y);
}



/*****************************************************************************************************/
/*****************************************************************************************************/

/*****************************************************************************************************/
/*****************************************************************************************************/

struct BranchState {
    vec<real> weights;
    std::unordered_set<int> keep_track;
    vec<std::pair<int, vec<int>>> stack;

    BranchState() = default;
    BranchState(NucSpace const &space) :
            weights(vmap(space.nucs, [](auto x) { return real(x.size()); })) {}

    auto const & top() const { return back(stack); }

    int distance_lower_bound(NucSpace const &space) {
        int count = 0;
        for (auto i : range(len(space.nucs))) {
            count += !(space.nucs[i].in((*space.ref)[i]));
        }
        return count;
    }

    int create_list(NucSpace &space, IV x, int i) {
        int min_dist = len(space.nucs);
        auto cur_dist = distance_lower_bound(space);

        vec<std::pair<int, int>> costs;

        for (int j = 0; j < 4; ++j) {
            if (x.in(j)) {
                std::unique_ptr<NucSpace> space_ptr(space.cast_clone());
                (*space_ptr).force(i, j);
                auto stat = (*space_ptr).status();

                int dist = len(space.nucs);
                if (stat != SpaceStatus::SS_FAILED) dist = distance_lower_bound(*space_ptr);
                if (dist < min_dist) min_dist = dist;
                costs.emplace_back(j, dist);
            }
        }

        std::sort(begin_of(costs), end_of(costs), [](auto const &a, auto const &b) {
            return a.second < b.second;
        });

        auto order = vmap(costs, [](auto const &x) { return x.first; });

        /* if nothing is workable, just give all values */
        if (!len(order)) {
            for (int j = 0; j < 4; ++j) {
                if (x.in(j)) order.emplace_back(j);
            }
        }

        /* update "stack" */
        keep_track.emplace(i);
        stack.emplace_back(std::make_pair(i, std::move(order)));

        /* update weight */
        weights[i] = 0.5 * weights[i] + (1 - 0.5) * (min_dist - cur_dist);
        // BEEP("just made", top().second);
        return top().second[0];
    }

    int choose_val(NucSpace &space, IV x, int i) {
        /* return first untried if already enumerated */
        if (keep_track.count(i)) {
            prune_back(i);
            // BEEP(top().second, x);
            for (auto j : top().second) if (x.in(j)) return j;
            NUPACK_BUG("should never get here!", top(), x);
        }

        return create_list(space, x, i);
    }

    void prune_back(int x) {
        while (top().first != x) {
            int i = top().first;
            // BEEP("popping: ", i);
            keep_track.erase(i);
            stack.pop_back();
        }
    }

    int last() const { return top().first; }

    auto size() const { return len(stack); }
};


class HeuristicBrancher : public Brancher {
protected:
    IVArray x;
    BranchState * state;

    // choice definition
    class PosVal : public Choice {
    public:
        int pos; int val;

        PosVal(const HeuristicBrancher& b, int p, int v)
            : Choice(b,2), pos(p), val(v) {}

        virtual void archive(Archive& e) const {
            Choice::archive(e);
            e << pos << val;
        }
    };

public:
    HeuristicBrancher(Home home, IVArray& x0, BranchState *_state)
        : Brancher(home), x(x0), state(_state) {}

    static void post(Home home, IVArray& x, BranchState * state) {
        (void) new (home) HeuristicBrancher(home, x, state);
    }

    virtual size_t dispose(Space& home) {
        (void) Brancher::dispose(home);
        return sizeof(*this);
    }

    HeuristicBrancher(Space& home, HeuristicBrancher& b)
        : Brancher(home,b) {
        x.update(home, b.x);
        state = b.state;
    }

    virtual Brancher* copy(Space& home) {
        return new (home) HeuristicBrancher(home, *this);
    }

    virtual bool status(const Space& home) const {
        for (int i=0; i<x.size(); i++)
            if (!x[i].assigned())
                return true;
        return false;
    }

    virtual Choice* choice(Space& home) {
        int var = -1;
        /* most recently added is fully determined */
        if (len(*state) && x[state->last()].size() != 1) var = state->last();
        else {
            vec<real> pos_vars;
            pos_vars.reserve(len(x));
            for (auto i: range(len(x))) {
                if (x[i].assigned()) pos_vars.emplace_back(0);
                else {
                    pos_vars.emplace_back(state->weights[i]);
                }
            }
            auto it = random_choice(pos_vars);
            var = it - begin_of(pos_vars);
        }

        // BEEP(var);

        auto val = state->choose_val(static_cast<NucSpace &>(home), x[var], var);
        return new PosVal(*this, var, val);

        GECODE_NEVER;
        return NULL;
    }

    virtual Choice* choice(const Space&, Archive& e) {
        int pos, val;
        e >> pos >> val;
        return new PosVal(*this, pos, val);
    }

    virtual ExecStatus commit(Space& home,
                                                        const Choice& c,
                                                        unsigned int a) {
        const PosVal& pv = static_cast<const PosVal&>(c);
        int pos=pv.pos, val=pv.val;
        if (a == 0)
            return me_failed(x[pos].eq(home,val)) ? ES_FAILED : ES_OK;
        else
            return me_failed(x[pos].nq(home,val)) ? ES_FAILED : ES_OK;
    }

    virtual void print(const Space& home, const Choice& c,
                                         unsigned int a,
                                         std::ostream& o) const {
        const PosVal& pv = static_cast<const PosVal&>(c);
        int pos=pv.pos, val=pv.val;
        if (a == 0)
            o << "x[" << pos << "] = " << val;
        else
            o << "x[" << pos << "] != " << val;
    }
};

void heuristic_branch(Home home, const IntVarArgs& x, BranchState * state) {
    if (home.failed()) return;
    IVArray y(home, x);
    HeuristicBrancher::post(home, y, state);
}
