/**
 * @brief Python-related C++ source code for rebind
 * @file Python.cc
 */
#include <rebind-python/API.h>
#include <rebind/Document.h>
#include <complex>
#include <any>
#include <iostream>

namespace rebind {

/******************************************************************************/

std::pair<Object, char const *> str(PyObject *o) {
    std::pair<Object, char const *> out({PyObject_Str(o), false}, nullptr);
    if (out.first) {
#       if PY_MAJOR_VERSION > 2
            out.second = PyUnicode_AsUTF8(out.first); // PyErr_Clear
#       else
            out.second = PyString_AsString(out.first);
#       endif
    }
    return out;
}

void print(PyObject *o) {
    auto p = str(o);
    if (p.second) std::cout << p.second << std::endl;
}

// Assuming a Python exception has been raised, fetch its string and put it in
// a C++ exception type. Does not clear the Python error status.
PythonError python_error(std::nullptr_t) noexcept {
    PyObject *type, *value, *traceback;
    PyErr_Fetch(&type, &value, &traceback);
    if (!type) return PythonError("Expected Python exception to be set");
    auto p = str(value);
    PyErr_Restore(type, value, traceback);
    return PythonError(p.second ? p.second : "Python error with failed str()");
}

/******************************************************************************/

/// type_index from PyBuffer format string (excludes constness)
std::type_info const & Buffer::format(std::string_view s) {
    auto it = std::find_if(Buffer::formats.begin(), Buffer::formats.end(),
        [&](auto const &p) {return p.first == s;});
    return it == Buffer::formats.end() ? typeid(void) : *it->second;
}

std::string_view Buffer::format(std::type_info const &t) {
    auto it = std::find_if(Buffer::formats.begin(), Buffer::formats.end(),
        [&](auto const &p) {return p.second == &t;});
    return it == Buffer::formats.end() ? std::string_view() : it->first;
}

std::size_t Buffer::itemsize(std::type_info const &t) {
    auto it = std::find_if(scalars.begin(), scalars.end(),
        [&](auto const &p) {return std::get<1>(p) == t;});
    return it == scalars.end() ? 0u : std::get<2>(*it) / CHAR_BIT;
}

/******************************************************************************/

std::string_view from_unicode(PyObject *o) {
    Py_ssize_t size;
#if PY_MAJOR_VERSION > 2
    char const *c = PyUnicode_AsUTF8AndSize(o, &size);
#else
    char *c;
    if (PyString_AsStringAndSize(o, &c, &size)) throw python_error();
#endif
    if (!c) throw python_error();
    return std::string_view(static_cast<char const *>(c), size);
}

std::string_view from_bytes(PyObject *o) {
    char *c;
    Py_ssize_t size;
    PyBytes_AsStringAndSize(+o, &c, &size);
    return std::string_view(c, size);
}

/******************************************************************************/

template <class T>
bool to_arithmetic(Object const &o, Variable &v) {
    DUMP("cast arithmetic in: ", v.type());
    if (PyFloat_Check(o)) return v = static_cast<T>(PyFloat_AsDouble(+o)), true;
    if (PyLong_Check(o)) return v = static_cast<T>(PyLong_AsLongLong(+o)), true;
    if (PyBool_Check(o)) return v = static_cast<T>(+o == Py_True), true;
    if (PyNumber_Check(+o)) { // This can be hit for e.g. numpy.int64
        if (std::is_integral_v<T>) {
            if (auto i = Object::from(PyNumber_Long(+o)))
                return v = static_cast<T>(PyLong_AsLongLong(+i)), true;
        } else {
            if (auto i = Object::from(PyNumber_Float(+o)))
               return v = static_cast<T>(PyFloat_AsDouble(+i)), true;
        }
    }
    DUMP("cast arithmetic out: ", v.type());
    return false;
}

/******************************************************************************/

bool object_response(Variable &v, TypeIndex t, Object o) {
    if (Debug) {
        auto repr = Object::from(PyObject_Repr(SubClass<PyTypeObject>{(+o)->ob_type}));
        DUMP("input object reference count", reference_count(o));
        DUMP("trying to convert object to ", t.name(), " ", from_unicode(+repr));
        DUMP(bool(cast_if<Variable>(o)));
    }

    if (auto p = cast_if<Variable>(o)) {
        DUMP("its a variable");
        Dispatch msg;
        v = p->request_variable(msg, t);
        return v.has_value();
    }

    if (t.matches<TypeIndex>()) {
        if (auto p = cast_if<TypeIndex>(o)) return v = *p, true;
        else return false;
    }

    if (t.equals<std::nullptr_t>()) {
        if (+o == Py_None) return v = nullptr, true;
    }

    if (t.matches<Function>()) {
        DUMP("requested function");
        if (+o == Py_None) v.emplace(Type<Function>());
        else if (auto p = cast_if<Function>(o)) v = *p;
        // general python function has no signature associated with it right now.
        // we could get them out via function.__annotations__ and process them into a tuple
        else v.emplace(Type<Function>())->emplace(PythonFunction({+o, true}, {Py_None, true}), {});
        return true;
    }

    if (t.equals<Sequence>()) {
        if (PyTuple_Check(o) || PyList_Check(o)) {
            DUMP("making a Sequence");
            Sequence *s = v.emplace(Type<Sequence>());
            s->reserve(PyObject_Length(o));
            map_iterable(o, [&](Object o) {s->emplace_back(std::move(o));});
            return true;
        } else return false;
    }

    if (t.equals<Real>())
        return to_arithmetic<Real>(o, v);

    if (t.equals<Integer>())
        return to_arithmetic<Integer>(o, v);

    if (t.equals<bool>()) {
        if ((+o)->ob_type == Py_None->ob_type) { // fix, doesnt work with Py_None...
            return v = false, true;
        } else return to_arithmetic<bool>(o, v);
    }

    if (t.equals<std::string_view>()) {
        if (PyUnicode_Check(+o)) return v.emplace(Type<std::string_view>(), from_unicode(+o)), true;
        if (PyBytes_Check(+o)) return v.emplace(Type<std::string_view>(), from_bytes(+o)), true;
        return false;
    }

    if (t.equals<std::string>()) {
        if (PyUnicode_Check(+o)) return v.emplace(Type<std::string>(), from_unicode(+o)), true;
        if (PyBytes_Check(+o)) return v.emplace(Type<std::string>(), from_bytes(+o)), true;
        return false;
    }

    if (t.equals<ArrayView>()) {
        if (PyObject_CheckBuffer(+o)) {
            // Read in the shape but ignore strides, suboffsets
            DUMP("cast buffer", reference_count(o));
            if (auto buff = Buffer(o, PyBUF_FULL_RO)) {
                DUMP("making data", reference_count(o));
                DUMP(Buffer::format(buff.view.format ? buff.view.format : "").name());
                DUMP("ndim", buff.view.ndim);
                DUMP((nullptr == buff.view.buf), bool(buff.view.readonly));
                for (auto i = 0; i != buff.view.ndim; ++i) DUMP(i, buff.view.shape[i], buff.view.strides[i]);
                DUMP("itemsize", buff.view.itemsize);
                ArrayLayout lay;
                lay.contents.reserve(buff.view.ndim);
                for (std::size_t i = 0; i != buff.view.ndim; ++i)
                    lay.contents.emplace_back(buff.view.shape[i], buff.view.strides[i] / buff.view.itemsize);
                DUMP("layout", lay, reference_count(o));
                DUMP("depth", lay.depth());
                ArrayData data{buff.view.buf, buff.view.format ? &Buffer::format(buff.view.format) : &typeid(void), !buff.view.readonly};
                return v.emplace(Type<ArrayView>(), std::move(data), std::move(lay)), true;
            } else throw python_error(type_error("C++: could not get buffer"));
        } else return false;
    }

    if (t.equals<std::complex<double>>()) {
        if (PyComplex_Check(+o)) return v.emplace(Type<std::complex<double>>(), PyComplex_RealAsDouble(+o), PyComplex_ImagAsDouble(+o)), true;
        return false;
    }

    DUMP("requested ", v.type(), t);
    return false;
}

/******************************************************************************/

std::string get_type_name(TypeIndex idx) noexcept {
    std::string out;
    auto it = type_names.find(idx);
    if (it == type_names.end() || it->second.empty()) out = idx.name();
    else out = it->second;
    out += QualifierSuffixes[static_cast<unsigned char>(idx.qualifier())];
    return out;
}

/******************************************************************************/

std::string wrong_type_message(WrongType const &e, std::string_view prefix) {
    std::ostringstream os;
    os << prefix << e.what() << " (#" << e.index << ", ";
    if (!e.source.empty())
        os << e.source << " \u2192 " << get_type_name(e.dest) << ", ";
    if (!e.indices.empty()) {
        auto it = e.indices.begin();
        os << "scopes=[" << *it;
        while (++it != e.indices.end()) os << ", " << *it;
        os << "], ";
    }
    if (e.expected != -1)
        os << "expected=" << e.expected << " received=" << e.received << ", ";
    std::string s = os.str();
    s.pop_back();
    s.back() = ')';
    return s;
}

/******************************************************************************/

Variable variable_reference_from_object(Object o) {
    if (auto p = cast_if<Function>(o)) return {Type<Function const &>(), *p};
    else if (auto p = cast_if<std::type_index>(o)) return {Type<std::type_index>(), *p};
    else if (auto p = cast_if<Variable>(o)) {
        DUMP("variable from object ", p, " ", p->data());
        DUMP("variable qualifier=", p->qualifier(), ", reference qualifier=", p->reference().qualifier());
        return p->reference();
    }
    else return std::move(o);
}

/******************************************************************************/

// Store the objects in args in pack
void args_from_python(Sequence &v, Object const &args) {
    v.reserve(v.size() + PyObject_Length(+args));
    map_iterable(args, [&v](Object o) {v.emplace_back(variable_reference_from_object(std::move(o)));});
}

/******************************************************************************/

}
