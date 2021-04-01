/**
 * @brief Python-related C++ API for rebind
 * @file PythonAPI.h
 */

#pragma once
#include "Object.h"

#include <rebind/Document.h>
#include <mutex>
#include <unordered_map>

namespace rebind {

extern std::unordered_map<TypeIndex, std::string> type_names;
extern Object TypeError, UnionType;
extern std::unordered_map<Object, Object> output_conversions, input_conversions, type_translations;
extern std::unordered_map<std::type_index, Object> python_types;

/******************************************************************************/

template <class ...Ts>
std::nullptr_t type_error(char const *s, Ts ...ts) {PyErr_Format(TypeError, s, ts...); return nullptr;}

struct Var : Variable {
    using Variable::Variable;
    Object ward = {};
};

template <>
struct Holder<Variable> : Holder<Var> {};

/******************************************************************************/

struct ArrayBuffer {
    std::vector<Py_ssize_t> shape_stride;
    std::size_t exports = 0;
    std::size_t n_elem;
    Object base;
    void *data;
    std::type_info const * type;
    bool mutate;

    ArrayBuffer() noexcept = default;
    ArrayBuffer(ArrayView const &a, Object const &b) : n_elem(a.layout.n_elem()),
        base(b), data(const_cast<void *>(a.data.pointer())), type(&a.data.type()), mutate(a.data.mutate()) {
        for (std::size_t i = 0; i != a.layout.depth(); ++i)
            shape_stride.emplace_back(a.layout.shape(i));
        auto const item = Buffer::itemsize(*type);
        for (std::size_t i = 0; i != a.layout.depth(); ++i)
            shape_stride.emplace_back(a.layout.stride(i) * item);
    }
};

/******************************************************************************/

Variable variable_reference_from_object(Object o);
void args_from_python(Sequence &s, Object const &pypack);
bool object_response(Variable &v, TypeIndex t, Object o);
std::string_view from_unicode(PyObject *o);

template <Qualifier Q>
struct Response<Object, Q> {
    bool operator()(Variable &v, TypeIndex t, Object o) const {
        DUMP("trying to get reference from qualified Object", Q, t);
        if (auto p = cast_if<Variable>(o)) {
            Dispatch msg;
            DUMP("requested qualified variable", t, p->type());
            v = p->reference().request_variable(msg, t);
            DUMP(p->type(), t, v.type());
        }
        return v.has_value();
    }
};

template <>
struct Response<Object, Value> {
    bool operator()(Variable &v, TypeIndex t, Object o) const {
        DUMP("trying to get reference from unqualified Object", t);
        if (!o) return false;
        DUMP("ref1", reference_count(o));
        Object type = Object(reinterpret_cast<PyObject *>((+o)->ob_type), true);
        if (auto p = input_conversions.find(type); p != input_conversions.end()) {
            Object guard(+o, false); // PyObject_CallFunctionObjArgs increments reference
            o = Object::from(PyObject_CallFunctionObjArgs(+p->second, +o, nullptr));
            type = Object(reinterpret_cast<PyObject *>((+o)->ob_type), true);
        }
        DUMP("ref2", reference_count(o));
        bool ok = object_response(v, t, std::move(o));
        DUMP("got response from object", ok);
        if (!ok) { // put diagnostic for the source type
            auto o = Object::from(PyObject_Repr(+type));
            DUMP("setting object error description", from_unicode(o));
            v = {Type<std::string>(), from_unicode(o)};
        }
        return ok;
    }
};

/******************************************************************************/

inline bool set_tuple_item(PyObject *t, Py_ssize_t i, PyObject *x) {
    if (!x) return false;
    incref(x);
    PyTuple_SET_ITEM(t, i, x);
    return true;
}

template <class V, class F>
Object map_as_tuple(V &&v, F &&f) noexcept {
    auto out = Object::from(PyTuple_New(std::size(v)));
    auto it = std::begin(v);
    using T = std::conditional_t<std::is_rvalue_reference_v<V>,
        decltype(*it), decltype(std::move(*it))>;
    for (Py_ssize_t i = 0u; i != std::size(v); ++i, ++it)
        if (!set_tuple_item(out, i, f(static_cast<T>(*it)))) return {};
    return out;
}

template <class ...Ts>
Object args_as_tuple(Ts &&...ts) {
    auto out = Object::from(PyTuple_New(sizeof...(Ts)));
    Py_ssize_t i = 0;
    auto go = [&](Object const &x) {return set_tuple_item(out, i++, x);};
    return (go(ts) && ...) ? out : Object();
}

Object variable_cast(Variable &&v, Object const &t={});

inline Object args_to_python(Sequence &&s, Object const &sig={}) {
    if (sig && !PyTuple_Check(+sig))
        throw python_error(type_error("expected tuple but got %R", (+sig)->ob_type));
    std::size_t len = sig ? PyTuple_GET_SIZE(+sig) : 0;
    auto const n = s.size();
    auto out = Object::from(PyTuple_New(n));
    Py_ssize_t i = 0u;
    for (auto &v : s) {
        if (i < len) {
            PyObject *t = PyTuple_GET_ITEM(+sig, i);
            throw python_error(type_error("conversion to python signature not implemented yet"));
        } else {
            // special case: if given an rvalue reference, make it into a value
            Variable &&var = v.qualifier() == Rvalue ? v.copy() : std::move(v);
            if (!set_tuple_item(out, i, variable_cast(std::move(var)))) return {};
        }
        ++i;
    }
    return out;
}

/******************************************************************************/

template <class F>
void map_iterable(Object iterable, F &&f) {
    auto iter = Object::from(PyObject_GetIter(+iterable));
    while (true) {
        if (auto it = Object(PyIter_Next(+iter), false)) f(std::move(it));
        else return;
    }
}

/******************************************************************************/

/// RAII release of Python GIL
struct PythonFrame final : Frame {
    std::mutex mutex;
    PyThreadState *state = nullptr;
    bool no_gil;

    PythonFrame(bool no_gil) : no_gil(no_gil) {}

    void enter() override {
        DUMP("running with nogil=", no_gil);
        if (no_gil && !state) state = PyEval_SaveThread(); // release GIL
    }

    std::shared_ptr<Frame> operator()(std::shared_ptr<Frame> &&t) override {
        DUMP("suspended Python ", bool(t));
        if (no_gil || state) return std::move(t); // return this
        else return std::make_shared<PythonFrame>(no_gil); // return a new frame
    }

    // acquire GIL; lock mutex to prevent multiple threads trying to get the thread going
    void acquire() noexcept {if (state) {mutex.lock(); PyEval_RestoreThread(state);}}
    // release GIL; unlock mutex
    void release() noexcept {if (state) {state = PyEval_SaveThread(); mutex.unlock();}}

    ~PythonFrame() {if (state) PyEval_RestoreThread(state);}
};

/******************************************************************************/

/// RAII reacquisition of Python GIL
struct ActivePython {
    PythonFrame &lock;

    ActivePython(PythonFrame &u) : lock(u) {lock.acquire();}
    ~ActivePython() {lock.release();}
};

/******************************************************************************/

struct PythonFunction {
    Object function, signature;

    PythonFunction(Object f, Object s={}) : function(std::move(f)), signature(std::move(s)) {
        if (+signature == Py_None) signature = Object();
        if (!function)
            throw python_error(type_error("cannot convert null object to Function"));
        if (!PyCallable_Check(+function))
            throw python_error(type_error("expected callable type but got %R", (+function)->ob_type));
        if (+signature && !PyTuple_Check(+signature))
            throw python_error(type_error("expected tuple or None but got %R", (+signature)->ob_type));
    }

    /// Run C++ functor; logs non-ClientError and rethrows all exceptions
    Variable operator()(Caller c, Sequence args) const {
        DUMP("calling python function");
        auto p = c.target<PythonFrame>();
        if (!p) throw DispatchError("Python context is expired or invalid");
        ActivePython lk(*p);
        Object o = args_to_python(std::move(args), signature);
        if (!o) throw python_error();
        return Variable(Object::from(PyObject_CallObject(function, o)));
    }
};

/******************************************************************************/

std::string get_type_name(TypeIndex idx) noexcept;

std::string wrong_type_message(WrongType const &e, std::string_view={});

namespace runtime {
    char const * unknown_exception_description() noexcept;
}

/******************************************************************************/

template <class F>
PyObject *raw_object(F &&f) noexcept {
    try {
        Object o = static_cast<F &&>(f)();
        xincref(+o);
        return +o;
    } catch (PythonError const &) {
        return nullptr;
    } catch (std::bad_alloc const &e) {
        PyErr_SetString(PyExc_MemoryError, "C++: out of memory (std::bad_alloc)");
    } catch (WrongNumber const &e) {
        unsigned int n0 = e.expected, n = e.received;
        PyErr_Format(TypeError, "C++: wrong number of arguments (expected %u, got %u)", n0, n);
    } catch (WrongType const &e) {
        try {PyErr_SetString(TypeError, wrong_type_message(e, "C++: ").c_str());}
        catch(...) {PyErr_SetString(TypeError, e.what());}
    } catch (std::exception const &e) {
        if (!PyErr_Occurred())
            PyErr_Format(PyExc_RuntimeError, "C++: %s", e.what());
    } catch (...) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, runtime::unknown_exception_description());
    }
    return nullptr;
}

/******************************************************************************/

template <class T>
PyObject * copy_from(PyObject *self, PyObject *other) noexcept {
    return raw_object([=] {
        cast_object<T>(self) = cast_object<T>(other); // not notexcept
        return Object(Py_None, true);
    });
}

/******************************************************************************/

}

