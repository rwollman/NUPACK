#pragma once
#include "Common.h"
#include "Error.h"

namespace rebind {

/******************************************************************************/

struct Dispatch;
struct VariableData;

using Storage = std::aligned_storage_t<4 * sizeof(void*), std::alignment_of_v<void *>>;

template <class T, class=void>
struct UseStack : std::integral_constant<bool, (sizeof(T) <= sizeof(Storage))
    && (std::alignment_of_v<Storage> % std::alignment_of_v<T> == 0)
    && std::is_nothrow_move_constructible_v<T>> {
    static_assert(std::is_same_v<T, std::decay_t<T>>);
};

/******************************************************************************/

enum class ActionType : std::uint_fast8_t {destroy, copy, move, response, assign};
using ActionFunction = void(*)(ActionType, void *, VariableData *);

struct RequestData {
    TypeIndex type;
    Dispatch *msg;
    Qualifier source;
};

static_assert(UseStack<RequestData>::value);
static_assert(std::is_trivially_destructible_v<RequestData>);

// destroy: delete the value stored at (void *)
// copy: copy value from (void *) into empty (Variable *)
// move: move value from (void *) into empty (Variable *)
// assign: assign existing value at (void *) = existing value in (Variable *)
// response: convert existing value in (void *) to new value in (Variable *)
//           using RequestData pre-stored in Variable->storage

/******************************************************************************/

// static_assert(sizeof(std::type_info const *) == sizeof(TypeIndex));
// static_assert(alignof(std::type_info const *) == alignof(TypeIndex));

struct VariableData {
    Storage buff; //< Buffer holding either pointer to the object, or the object itself
    ActionFunction act; //< Action<T>::apply of the held object, or NULL
    TypeIndex idx; //< type and qualifier of the held object, or NULL
    bool stack; //< Whether the held type (non-reference) can fit in the buffer

    constexpr VariableData() noexcept : buff(), act(nullptr), stack(false) {}

    VariableData(TypeIndex i, ActionFunction a, bool s) noexcept : buff(), act(a), idx(i), stack(s) {}

    void reset_data() noexcept {
        if (!act) return;
        buff = Storage();
        idx = {};
        act = nullptr;
        stack = false;
    }

    /// Return a pointer to the held object if it exists
    void * pointer() const noexcept {
        if (!act) return nullptr;
        else if (idx.qualifier() == Value && stack) return &const_cast<Storage &>(buff);
        else return reinterpret_cast<void * const &>(buff);
    }

    /// Return a pointer to the held object if it exists and is being managed
    void * handle() const noexcept {
        if (!act || idx.qualifier() != Value) return nullptr;
        else if (stack) return &const_cast<Storage &>(buff);
        else return reinterpret_cast<void * const &>(buff);
    }

    template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
    std::remove_reference_t<T> *target_pointer(Type<T> t, Qualifier q) const noexcept {
        // Qualifier is assumed not to be V
        return (idx.matches<T>()) && (
            (std::is_const_v<std::remove_reference_t<T>>)
            || (std::is_rvalue_reference_v<T> && q == Rvalue)
            || (std::is_lvalue_reference_v<T> && (q == Lvalue))
        ) ? static_cast<std::remove_reference_t<T> *>(pointer()) : nullptr;
    }
};

/******************************************************************************/

}
