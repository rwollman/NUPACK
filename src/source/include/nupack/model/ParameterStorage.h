#pragma once

namespace nupack {

/******************************************************************************************/

template <std::size_t I, std::size_t ...Is>
struct ParameterArray {
    static constexpr std::size_t begin = I;
    static constexpr std::size_t end = begin + (1 * ... * Is);
    static constexpr std::size_t ndim = sizeof...(Is);
    constexpr std::size_t size() const {return end - begin;}
    constexpr std::size_t back() const {return size() - 1;}

    template <class ...Js>
    static std::size_t index(Js ...js) {
        static_assert(sizeof...(Js) == sizeof...(Is));
        std::size_t out = begin, stride = 1;
        ((out += js * stride, stride *= Is), ...); // column major
        return out;
    }

    template <class T, std::size_t ...Js>
    static std::size_t array_index(T const &t, std::index_sequence<Js...>) {return index(t[Js]...);}

    template <class T>
    static std::size_t array_index(T const &t) {return array_index(t, std::make_index_sequence<ndim>());}
};

template <std::size_t I, std::size_t ...Is>
ParameterArray<I, (Is ? 4 : 4)...> parameter_grid(std::index_sequence<Is...>); //undefined

template <std::size_t I, std::size_t N>
using ParameterGrid = decltype(parameter_grid<I>(std::make_index_sequence<N>()));

/******************************************************************************************/

static constexpr ParameterGrid <0, 8> interior_2_2;
static constexpr ParameterGrid <decltype(interior_2_2)::end, 7> interior_1_2;
static constexpr ParameterGrid <decltype(interior_1_2)::end, 6> interior_1_1;
static constexpr ParameterGrid <decltype(interior_1_1)::end, 4> interior_mismatch;
static constexpr ParameterGrid <decltype(interior_mismatch)::end, 4> terminal_mismatch;
static constexpr ParameterGrid <decltype(terminal_mismatch)::end, 4> stack;
static constexpr ParameterGrid <decltype(stack)::end, 4> coaxial_stack;
static constexpr ParameterGrid <decltype(coaxial_stack)::end, 6> hairpin_tetra;
static constexpr ParameterGrid <decltype(hairpin_tetra)::end, 5> hairpin_tri;
static constexpr ParameterGrid <decltype(hairpin_tri)::end, 4> hairpin_mismatch;
static constexpr ParameterGrid <decltype(hairpin_mismatch)::end, 3> dangle5;
static constexpr ParameterGrid <decltype(dangle5)::end, 3> dangle3;
static constexpr ParameterGrid <decltype(dangle3)::end, 2> terminal_penalty;
static constexpr ParameterArray<decltype(terminal_penalty)::end, 30> interior_size;
static constexpr ParameterArray<decltype(interior_size)::end, 30> bulge_size;
static constexpr ParameterArray<decltype(bulge_size)::end, 30> hairpin_size;
static constexpr ParameterArray<decltype(hairpin_size)::end, 5> ninio;
static constexpr ParameterArray<decltype(ninio)::end> multi_base;
static constexpr ParameterArray<decltype(multi_base)::end> multi_init;
static constexpr ParameterArray<decltype(multi_init)::end> multi_pair;
static constexpr ParameterArray<decltype(multi_pair)::end> log_loop_penalty;
static constexpr ParameterArray<decltype(log_loop_penalty)::end> join_penalty;

/******************************************************************************************/

}
