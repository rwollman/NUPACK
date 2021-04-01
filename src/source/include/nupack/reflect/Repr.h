#pragma once
#include "../algorithms/Operators.h"
#include "../algorithms/Traits.h"
#include "../standard/Ptr.h"
#include <atomic>

namespace nupack {

NUPACK_DETECT(has_repr, decltype(declref<T const>().repr()));
/// repr is a unique representation
template <class T, NUPACK_IF(has_repr<T>)>
auto repr_of(T const &t) {return t.repr();}
template <class T, NUPACK_IF(is_empty<T> && !has_repr<T>)>
constexpr auto repr_of(T const &t) {return string();}

/******************************************************************************************/
extern std::atomic<std::size_t> unique_repr_count;
inline auto unique_count() {return std::to_string(unique_repr_count++);}

template <class T, NUPACK_IF(!is_empty<T> && !has_repr<T>)>
auto repr_of(T const &t) {return unique_count();}

/// for non-introspected classes we just keep a counter for the objects
class Unique_Repr : public TotallyOrdered {
    Optional<std::type_index> type;
    string data;

public:

    NUPACK_REFLECT(Unique_Repr, type, data);

    friend std::ostream & operator<<(std::ostream &os, Unique_Repr const &) {return os << "Unique_Repr()";}

    Unique_Repr() = default;

    template <class T>
    Unique_Repr(T const &t) : type(typeid(T)), data(repr_of(t)) {}

    bool operator<(Unique_Repr const &o) const {
        if (!type || !o.type) return false;
        return std::tie(type, data) < std::tie(o.type, o.data);
    }

    bool operator==(Unique_Repr const &o) const {
        if (!type || !o.type) return true; // any matches anything
        return std::tie(*type, data) == std::tie(*o.type, o.data);
    }
};


}
