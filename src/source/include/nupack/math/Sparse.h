#pragma once
#include "../types/Matrix.h"

namespace nupack {

/******************************************************************************************/

template <class T>
real sparsity(T const &t) {return 1 - std::count(begin_of(t), end_of(t), 0) / real(len(t));}

/******************************************************************************************/

/// Make a matrix from a list of tuples where it is nonzero (x, y, value)
template <class Mat, class Tuples, NUPACK_IF(la::is_dense<Mat>)>
Mat matrix_from_tuples(usize n_rows, usize n_cols, Tuples const &tups) {
    Mat ret(n_rows, n_cols); la::fill_zero(ret);
    for (auto const &t : tups) ret(first_of(t), second_of(t)) += third_of(t);
    return ret;
};

/******************************************************************************************/

template <class Mat, class Tuples, NUPACK_IF(la::is_sparse<Mat>)>
Mat matrix_from_tuples(usize const n_rows, usize const n_cols, Tuples const &tups) {
    la::vec values(len(tups));
    la::umat locs(2, len(tups));
    for (auto i : indices(tups)) {
        locs(0, i) = first_of(tups[i]);
        locs(1, i) = second_of(tups[i]);
        values(i) = third_of(tups[i]);
    }
    return {locs, values, n_rows, n_cols};
}

/******************************************************************************************/

template <class T>
struct SparsePairs {
    Col<T> values, diag;
    Col<std::uint32_t> rows, cols;

    NUPACK_REFLECT(SparsePairs, values, diag, rows, cols);
};

/// Take square matrix of probabilities
/// Return diagonal entries, off-diagonal entries, rows, columns
/// if n is 0, return all entries, else return the n highest non-diagonal entries in each row
// complexity is N^2 log(n) + N n log(n N) \approx N^2 log(n) + N n log(N)
template <class T>
auto sparse_pair_matrix(Mat<T> const &m, std::size_t const row_size=0, T threshold=0) {
    NUPACK_REQUIRE(m.n_rows, ==, m.n_cols);
    NUPACK_ASSERT(m.is_symmetric(), "sparse_pair_matrix(): pair matrix should be symmetric", m);
    SparsePairs<T> o;
    o.diag = m.diag();

    if (m.n_rows <= 1) return o;

    auto const filter = [threshold](T const &t) {return t > threshold;};
    auto const nnz = (count_if(m, filter) - count_if(o.diag, filter)) / 2;
    bool const simple = (row_size == 0) || (row_size >= m.n_rows/2) || (nnz <= row_size * m.n_rows);

    if (simple) {
        for_each(std::tie(o.values, o.rows, o.cols), [nnz](auto &x) {x.set_size(nnz);});
        la::uword p = 0;
        for (auto j : range(m.n_rows)) for (auto i : range(0, j)) {
            if (filter(m.at(i, j))) {
                o.rows(p) = i;
                o.cols(p) = j;
                o.values(p) = m.at(i, j);
                ++p;
            }
        }
        NUPACK_REQUIRE(p, ==, nnz, "Sparse pair matrix failure", count_if(m, filter), count_if(o.diag, filter));
    } else {
        std::vector<la::uword> idx;
        std::vector<std::pair<la::uword, la::uword>> v;
        v.reserve(m.n_rows * row_size);

        for (auto const j : range(m.n_rows)) {
            // all possible column indices
            idx.assign(range(m.n_rows).begin(), range(m.n_rows).end());
            // partially sort to get the n biggest column indices
            std::nth_element(idx.begin(), idx.begin() + row_size, idx.end(), [j, c=m.col(j)](auto i1, auto i2) {
                if (i1 == j) return false; // don't include on-diagonal elements
                if (i2 == j) return true;
                return c(i1) > c(i2);
            });
            // only include nonzero elements
            for (auto it = idx.begin(); it != idx.begin() + row_size; ++it)
                if (filter(m.at(*it, j))) v.emplace_back(std::minmax(*it, j));
        }
        v = unique_sorted(std::move(v));

        for_each(std::tie(o.values, o.rows, o.cols), [s=v.size()](auto &x) {x.set_size(s);});
        zip(v, o.rows, o.cols, o.values, [&m](auto const &p, auto &r, auto &c, auto &v) {
            r = p.first;
            c = p.second;
            v = m.at(p.first, p.second);
        });
    }

    return o;
}


}
