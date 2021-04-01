#pragma once

#include "../types/Sequence.h"
#include "../reflect/Print.h"
#include "../iteration/Search.h"
#include "../iteration/Transform.h"

namespace nupack {

/****************************************************************************************/

/// EdgeSet describes a node on graph
struct EdgeSet : ConstIndexable<EdgeSet>, TotallyOrdered {
    using size_type = iseq;
    using value_type = EdgeList::value_type;
    using const_iterator = EdgeList::const_iterator;

    /************************************************************************************/

    NUPACK_REFLECT(EdgeSet, index, parent, edges, parent_loc);

    auto & iter() const {return edges;}

    /// All neighbors, ordered (with rotation symmetry)
    EdgeList edges;
    /// This EdgeSet's index, and the index of its parent
    Edge index, parent;
    /// Location of parent within edges
    size_type parent_loc;

    /************************************************************************************/

    EdgeSet() = default;

    EdgeSet(Edge i, Edge p) : index(i), parent(p), parent_loc(0) {edges.emplace_back(p);};
    EdgeSet(Edge i, EdgeList e, Edge p) : edges(std::move(e)), index(i), parent(p) {update_parent_loc();};

    /************************************************************************************/

    /// Does this EdgeSet have a root?
    bool is_root() const;
    /// Check internal index consistency
    bool check() const;

    /************************************************************************************/

    /// Change parent to a new neighbor, thus changing its location
    void set_parent(Edge p) {parent = p; update_parent_loc();}
    /// Replace parent index with a new one, location not changing
    void replace_parent(Edge p) {parent = p; edges[parent_loc] = p;}
    /// Transfer to new index
    template <class W> void transfer(Edge, W const &);
    /// Replace neighbor
    void replace(Edge from, Edge to);
    /// Switch the parent
    Edge flip(Edge new_parent);
    /// Rotate EdgeSet and update root location accordingly
    void rotate(size_type shift);

    void append(Edge e) {edges.emplace_back(e);}

    /************************************************************************************/

    /// Make sure parent_loc points to parent
    void update_parent_loc();

    /************************************************************************************/

    /// Merge this loop with its child
    template <bool Check=true, class W> std::pair<Edge, Edge> merge(EdgeSet const &k, W const &);
    /// Merge specialized for two exterior loops (results in dissociated complex)
    template <class W> std::pair<Edge, Edge> dissociate(EdgeSet &k, int, int, W const &);
    /// Join two exterior loops (results in associated complex)
    template <class W> void associate(EdgeSet &k, size_type, size_type, int, int, W const &);
    /// Make a new edge, using its index, starting and ending indices in this set
    template <class W> EdgeSet split(Edge new_self, size_type s1, size_type s2, W const &);
    /// Get indices of parent in kid loop, kid in parent loop
    static std::pair<size_type, size_type> get_locs(EdgeSet const &, EdgeSet const &);

    /************************************************************************************/

    bool operator==(EdgeSet const &e) const {
        return std::tie(index, parent, edges) == std::tie(e.index, e.parent, e.edges);
    }

    bool operator<(EdgeSet const &e) const {
        return std::tie(index, parent, edges) < std::tie(e.index, e.parent, e.edges);
    }

    /************************************************************************************/
    /// Get edge const_iterator from edge value
    auto find_edge(Edge e) const {
        auto out = find(edges, e);
        NUPACK_DASSERT(out != end_of(edges), "Could not find edge (const)");
        return out;
    }
    /// Get edge iterator from edge value
    auto find_edge(Edge e) {
        auto out = find(edges, e);
        NUPACK_ASSERT(out != end_of(edges), "Could not find edge");
        return out;
    }
    /// Get edge index from edge value
    auto find_edge_index(Edge e) const {return find_edge(e) - begin_of(edges);}
};

/****************************************************************************************/

inline bool EdgeSet::is_root() const {
    NUPACK_DREQUIRE(edges.at(parent_loc), ==, parent);
    return edges[parent_loc] == Ether;
}

/****************************************************************************************/

inline bool EdgeSet::check() const {
    if (!contains(edges, parent)) return false;
    usize sum = 0;
    for (auto i : edges) for (auto j : edges) sum += (j == i);
    return sum == len(edges);
}

/****************************************************************************************/

template <class W> void EdgeSet::transfer(Edge i, W const &w) {
    for (auto const &e : edges) if (e != Ether && e != parent)
        w(e).replace_parent(i);
    if (parent != Ether) w(parent).replace(index, i);
    index = i;
}

/****************************************************************************************/

/// Find the index of the first in the second, second in the first
inline auto EdgeSet::get_locs(EdgeSet const &e1, EdgeSet const &e2) -> std::pair<size_type, size_type>{
    if (!Debug) return {e2.find_edge_index(e1.index), e1.find_edge_index(e2.index)};
    try {return {e2.find_edge_index(e1.index), e1.find_edge_index(e2.index)};}
    catch(...) {NUPACK_BUG("Failure in get_locs", e1, e1.index, e2, e2.index);}
}

/****************************************************************************************/

inline void EdgeSet::replace(Edge from, Edge to) {*find_edge(from) = to;}

/****************************************************************************************/

inline Edge EdgeSet::flip(Edge new_parent) {
    auto out = edges[parent_loc];
    parent_loc = find_edge(new_parent) - begin_of(edges);
    return out;
}

/****************************************************************************************/

inline void EdgeSet::rotate(size_type shift) {
    NUPACK_DREQUIRE(shift, <, len(edges));
    if (!shift) return;
    std::rotate(begin_of(edges), begin_of(edges) + shift, end_of(edges));
    parent_loc += len(edges) - shift;
    if (parent_loc >= len(edges)) parent_loc -= len(edges);
}

/****************************************************************************************/

inline void EdgeSet::update_parent_loc() {
    parent_loc = find_edge(parent) - begin_of(edges);
}

/****************************************************************************************/

template <bool Check, class W>
std::pair<Edge, Edge> EdgeSet::merge(EdgeSet const &k, W const &w) {
    NUPACK_DREQUIRE(k.parent_loc, ==, k.find_edge_index(index));

    auto kp = find_edge(k.index);
    auto pk = k.begin() + k.parent_loc;

    for (auto e : k.edges) if (e != index && e != Ether) w(e).replace_parent(index);

    std::pair<Edge, Edge> out = {pk - k.begin(), kp - begin()};

    edges = catted<EdgeList>(kp + 1, end(), begin(), kp, pk + 1, k.end(), k.begin(), pk);

    update_parent_loc();

    if (Check) for (auto e : edges) if (e != parent && e != Ether)
        NUPACK_DREQUIRE(w(e).parent, ==, index, e);

    return out;
}

/****************************************************************************************/

template <class W>
std::pair<Edge, Edge> EdgeSet::dissociate(EdgeSet &k, int nick, int k_nick, W const &w) {
    auto kp = find_edge(k.index), pk = k.find_edge(index);
    auto out = std::make_pair(pk - k.begin(), kp - begin());

    EdgeList p_edges = {Ether}, k_edges = {Ether};
    circular_cat(p_edges, edges, begin() + nick + 1, kp);
    circular_cat(p_edges, k.edges, pk + 1, k.begin() + k_nick);

    circular_cat(k_edges, k.edges, k.begin() + k_nick + 1, pk);
    circular_cat(k_edges, edges, kp + 1, begin() + nick);

    edges = std::move(p_edges); k.edges = std::move(k_edges);


    if (contains(edges, parent)) {
        k.set_parent(Ether);
        update_parent_loc();
    } else {
        k.set_parent(parent);
        if (parent != Ether) w(parent).replace(index, k.index);
        set_parent(Ether);
    }

    for (auto const &e : k.edges) if (e != Ether && e != k.parent) w(e).replace_parent(k.index);
    for (auto const &e : edges) if (e != Ether && e != parent) w(e).replace_parent(index);

    return out;
}

/****************************************************************************************/

template <class W>
void EdgeSet::associate(EdgeSet &k, size_type s, size_type ks, int nick, int k_nick, W const &w) {
    {NUPACK_DREQUIRE(this, !=, &k); NUPACK_DASSERT(is_root());}

    EdgeList p_edges = {k.index};
    circular_cat(p_edges, edges, begin() + s + 1, begin() + nick);
    p_edges.emplace_back(Ether);
    circular_cat(p_edges, k.edges, k.begin() + k_nick + 1, k.begin() + ks + 1);

    EdgeList k_edges = {index};
    circular_cat(k_edges, k.edges, k.begin() + ks + 1, k.begin() + k_nick);
    k_edges.emplace_back(Ether);
    circular_cat(k_edges, edges, begin() + nick + 1, begin() + s + 1);

    edges = std::move(p_edges); k.edges = std::move(k_edges);

    if (contains(k.edges, k.parent)) {
        k.update_parent_loc();
        set_parent(k.index);
    } else {
        set_parent(k.parent);
        if (parent != Ether) w(parent).replace(k.index, index);
        k.set_parent(index);
    }
    for (auto const &e : k.edges) if (e != Ether && e != k.parent) w(e).replace_parent(k.index);
    for (auto const &e : edges) if (e != Ether && e != parent) w(e).replace_parent(index);
}

/****************************************************************************************/

template <class W> EdgeSet EdgeSet::split(Edge new_self, size_type s1, size_type s2, W const &w) {
    NUPACK_DREQUIRE(s1, <=, s2);
    NUPACK_DREQUIRE(s2, <=, len(edges));
    EdgeSet out(new_self, index);

    // Get edges of new edge set
    extend(out.edges, begin() + s1, begin() + s2);
    // Update children of new edge set to their new parent
    for (auto e : range(out.begin() + 1, out.end())) {
        if (*e != Ether && *e != parent) w(*e).replace_parent(out.index);
    }

    auto it = edges.erase(begin_of(edges) + s1, begin_of(edges) + s2);
    edges.insert(it, out.index);

    if (contains(edges, parent)) update_parent_loc();
    else {
        out.set_parent(parent);
        if (parent != Ether) w(parent).replace(index, out.index);
        set_parent(out.index);
    }

    return out;
}

/****************************************************************************************/

}
