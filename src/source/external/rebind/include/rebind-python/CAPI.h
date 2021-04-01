#pragma once
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wregister"
#include <Python.h>
#pragma GCC diagnostic pop

#include <functional>
#include <rebind/Type.h>
#include <rebind/Common.h>
#include <rebind/Error.h>

namespace rebind {

/******************************************************************************/

inline std::size_t reference_count(PyObject *o) {return o ? Py_REFCNT(o) : 0u;}
inline void incref(PyObject *o) noexcept {Py_INCREF(o);}
inline void decref(PyObject *o) noexcept {Py_DECREF(o);}
inline void xincref(PyObject *o) noexcept {Py_XINCREF(o);}
inline void xdecref(PyObject *o) noexcept {Py_XDECREF(o);}
inline PyObject * not_none(PyObject *o) {return o == Py_None ? nullptr : o;}

/******************************************************************************/

void print(PyObject *o);

using Version = std::tuple<unsigned, unsigned, unsigned>;
static constexpr Version PythonVersion{PY_MAJOR_VERSION, PY_MINOR_VERSION, PY_MICRO_VERSION};

/******************************************************************************/

template <class T>
struct SubClass {
    T *ptr;
    operator PyObject *() const {return reinterpret_cast<PyObject *>(ptr);}
    operator T *() const {return ptr;}
};

template <class T>
struct Holder {
    static inline PyTypeObject type;
    PyObject_HEAD // 16 bytes for the ref count and the type object
    T value; // I think stack is OK because this object is only casted to anyway.
};

template <class T>
SubClass<PyTypeObject> type_object(Type<T> t={}) {return {&Holder<T>::type};}

/******************************************************************************/

template <class T>
T * cast_if(PyObject *o) {
    if (!PyObject_TypeCheck(o, type_object<T>())) return nullptr;
    return std::addressof(reinterpret_cast<Holder<T> *>(o)->value);
}

template <class T>
T & cast_object(PyObject *o) {
    if (!PyObject_TypeCheck(o, type_object<T>()))
        throw std::invalid_argument("Expected instance of rebind.TypeIndex");
    return reinterpret_cast<Holder<T> *>(o)->value;
}

/******************************************************************************/

template <class T>
PyObject *tp_new(PyTypeObject *subtype, PyObject *, PyObject *) noexcept {
    static_assert(noexcept(T{}), "Default constructor should be noexcept");
    PyObject *o = subtype->tp_alloc(subtype, 0); // 0 unused
    if (o) new (&cast_object<T>(o)) T; // Default construct the C++ type
    return o;
}

template <class T>
void tp_delete(PyObject *o) noexcept {
    reinterpret_cast<Holder<T> *>(o)->~Holder<T>();
    Py_TYPE(o)->tp_free(o);
}

/******************************************************************************/

class Buffer {
    static Vector<std::pair<std::string_view, std::type_info const *>> formats;
    bool valid;

public:
    Py_buffer view;

    Buffer(Buffer const &) = delete;
    Buffer(Buffer &&b) noexcept : view(b.view), valid(std::exchange(b.valid, false)) {}

    Buffer & operator=(Buffer const &) = delete;
    Buffer & operator=(Buffer &&b) noexcept {view = b.view; valid = std::exchange(b.valid, false); return *this;}

    explicit operator bool() const {return valid;}

    Buffer(PyObject *o, int flags) {
        DUMP("before buffer", reference_count(o));
        valid = PyObject_GetBuffer(o, &view, flags) == 0;
        if (valid) DUMP("after buffer", reference_count(o), view.obj == o);
    }

    static std::type_info const & format(std::string_view s);
    static std::string_view format(std::type_info const &t);
    static std::size_t itemsize(std::type_info const &t);
    // static Binary binary(Py_buffer *view, std::size_t len);
    // static Variable binary_view(Py_buffer *view, std::size_t len);

    ~Buffer() {
        if (valid) {
            PyObject *o = view.obj;
            DUMP("before release", reference_count(view.obj), " ", view.obj);
            PyBuffer_Release(&view);
            DUMP("after release", reference_count(o));
        }
    }
};

/******************************************************************************/

template <class T, class U>
bool compare(decltype(Py_EQ) op, T const &t, U const &u) {
    switch(op) {
        case(Py_EQ): return t == u;
        case(Py_NE): return t != u;
        case(Py_LT): return t <  u;
        case(Py_GT): return t >  u;
        case(Py_LE): return t <= u;
        case(Py_GE): return t >= u;
    }
    return false;
}

/******************************************************************************/

template <class T>
PyTypeObject type_definition(char const *name, char const *doc) {
    PyTypeObject o{PyVarObject_HEAD_INIT(NULL, 0)};
    o.tp_name = name;
    o.tp_basicsize = sizeof(Holder<T>);
    o.tp_dealloc = tp_delete<T>;
    o.tp_new = tp_new<T>;
    o.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    o.tp_doc = doc;
    return o;
}

}
