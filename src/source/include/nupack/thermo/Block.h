/**
 * @brief Generic Block CRTP object providing compile-time functionality
 *
 * @file Block.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Adapters.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

template <class Contents, class Names>
struct Record {
    Contents contents;
    using value_type = value_type_of<std::tuple_element_t<1, Contents>>;

    template <usize ...Is>
    auto members_(indices_t<Is...>) {return make_members(at<Is>(contents)...);}
    auto members() {return members_(indices_in(contents));}
    NUPACK_EXTEND_NAMES(Names, complete);
    bool complete() const {return at<tuple_size<Contents>-1>(contents);}
};

NUPACK_DEFINE_TEMPLATE(is_record, Record, class, class);

/******************************************************************************************/

// Block takes the type T and a template to fill with it
// for T & inherit from Base<T>
// for T inherit from Base<T> AND Base<T>::storage_type

template <class O, class M, usize ...Is>
O get_subsquare(M &&m, span s, indices_t<Is...>) {return{True(), at<Is>(m).subsquare(s)...};}


template <class T, class Base>
struct Block : Base::template storage_type<T>, Base, MemberComparable {
    using value_type = decay<T>;
    using base_type = typename Base::template storage_type<T>;

private:
    /// Indices helper function
    static constexpr auto indices() {return member_indices<base_type>();}

    template <usize ...Is>
    auto write(span i, span j, bool complete, indices_t<Is...>) const {
        static_assert(!is_ref<T>);
        return std::make_tuple(at<Is>(members_of(*this)).write(i, j, complete)..., complete);
    }

    template <class Ts, usize ...Is>
    void read(span i, span j, Ts const &ts, indices_t<Is...>) {
        static_assert(!is_ref<T>);
        NUPACK_UNPACK(at<Is>(members_of(*this)).read(i, j, at<Is>(ts)));
    }

    template <class B, usize ...Is>
    Block(False, indices_t<Is...>, B const &b) : base_type{at<Is>(members_of(b))...} {}

    template <class B, usize ...Is>
    Block(False, indices_t<Is...>, B &&b) : base_type{std::move(at<Is>(members_of(b)))...} {}

public:
    template <class U, class ...Ts, NUPACK_IF(!is_like<U, Block, True, False>)>
    Block(U &&u, Ts &&...ts) : base_type{Base::template storage<T>(fw<U>(u), fw<Ts>(ts)...)} {}

    template <class ...Ts>
    Block(True, Ts &&...ts) : base_type{fw<Ts>(ts)...} {}

    template <class T2, class B2>
    Block(Block<T2, B2> const &b) : Block(False(), indices(), b) {static_assert(!is_same<T, T2>, "");}

    template <class T2, class B2>
    Block(Block<T2, B2> &&b) : Block(False(), indices(), std::move(b)) {static_assert(!is_same<T, T2>, "");}

    template <class Record>
    void read(span i, span j, Record const &r) {read(i, j, r.contents, indices());}

    auto write(span i, span j, bool complete) const {
        using base = typename Base::template storage_type<value_type>;
        return Record<decltype(write(i, j, complete, indices())), base>{write(i, j, complete, indices())};
    }

    /// Get square submatrix for each member
    auto subsquare(span s) const {return get_subsquare<Block<value_type const &, Base>>(members_of(*this), s, indices());}
    auto subsquare(span s) {return get_subsquare<Block<value_type &, Base>>(members_of(*this), s, indices());}
    /// Return length of the matrix - I lazily use 1 not 0 because X doesn't have size
    auto size() const {return len(at<1>(members_of(*this)));}
    /// Return highest level result in this block
    auto result() const {return value_of(base_type::Q(0, size() - 1));}

    void copy_square(span i, span j) {
        for_each(members_of(*this), [=](auto &M) {M.copy_square(i, j);});
    }
};

NUPACK_DEFINE_TEMPLATE(is_block, Block, class, class);


/******************************************************************************************/

}}
