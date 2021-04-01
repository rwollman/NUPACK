#include <nupack/design/DesignComponents.h>

namespace nupack { namespace newdesign {

vec<real> ord_lin_lsq(vec<real> const &x, vec<real> const &y) {
    NUPACK_REQUIRE(len(x), ==, len(y));

    real_col col_y(y);

    real_mat X(len(x), 2, arma::fill::ones);
    real_col col_x(x);
    X.col(1) = col_x;
    auto Xt = X.t();

    real_col beta = solve(Xt * X, Xt * col_y);
    NUPACK_REQUIRE(len(beta), ==, 2);
    return {beta.begin(), beta.end()};
}

Timer & Timer::start() {
    _start = std::chrono::high_resolution_clock::now();
    _stop = _start;
    return *this;
}

real Timer::elapsed() const {
    return std::chrono::duration<real>(std::chrono::high_resolution_clock::now() - _start).count();
}

real Timer::stop() {
    if (!(_stop > _start)) _stop = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<real>(_stop - _start).count();
}

}}
