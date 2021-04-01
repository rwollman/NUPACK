#include <rebind-python/Cast.h>
#include <numeric>

namespace rebind {

bool is_subclass(PyTypeObject *o, PyTypeObject *t) {
    int x = PyObject_IsSubclass(reinterpret_cast<PyObject *>(o), reinterpret_cast<PyObject *>(t));
    return (x < 0) ? throw python_error() : x;
}

/******************************************************************************/

Object type_args(Object const &o) {
    auto out = Object::from(PyObject_GetAttrString(o, "__args__"));
    if (out && !PyTuple_Check(out))
        throw python_error(type_error("expected __args__ to be a tuple"));
    return out;
}

Object type_args(Object const &o, Py_ssize_t n) {
    auto out = type_args(o);
    if (out) {
        Py_ssize_t const m = PyTuple_GET_SIZE(+out);
        if (m != n) throw python_error(type_error("expected __args__ to be length %zd (got %zd)", n, m));
    }
    return out;
}

Object list_cast(Variable &&ref, Object const &o, Object const &root) {
    DUMP("Cast to list ", ref.type());
    if (auto args = type_args(o, 1)) {
        DUMP("is list ", PyList_Check(o));
        auto v = ref.cast<Sequence>();
        Object vt{PyTuple_GET_ITEM(+args, 0), true};
        auto list = Object::from(PyList_New(v.size()));
        for (Py_ssize_t i = 0; i != v.size(); ++i) {
            DUMP("list index ", i);
            Object item = python_cast(std::move(v[i]), vt, root);
            if (!item) return {};
            incref(+item);
            PyList_SET_ITEM(+list, i, +item);
        }
        return list;
    } else return {};
}

Object tuple_cast(Variable &&ref, Object const &o, Object const &root) {
    DUMP("Cast to tuple ", ref.type());
    if (auto args = type_args(o)) {
        Py_ssize_t const len = PyTuple_GET_SIZE(+args);
        auto v = ref.cast<Sequence>();
        if (len == 2 && PyTuple_GET_ITEM(+args, 1) == Py_Ellipsis) {
            Object vt{PyTuple_GET_ITEM(+args, 0), true};
            auto tup = Object::from(PyTuple_New(v.size()));
            for (Py_ssize_t i = 0; i != v.size(); ++i)
                if (!set_tuple_item(tup, i, python_cast(std::move(v[i]), vt, root))) return {};
            return tup;
        } else if (len == v.size()) {
            auto tup = Object::from(PyTuple_New(len));
            for (Py_ssize_t i = 0; i != len; ++i)
                if (!set_tuple_item(tup, i, python_cast(std::move(v[i]), {PyTuple_GET_ITEM(+args, i), true}, root))) return {};
            return tup;
        }
    }
    return {};
}

Object dict_cast(Variable &&ref, Object const &o, Object const &root) {
    DUMP("Cast to dict ", ref.type());
    if (auto args = type_args(o, 2)) {
        Object key{PyTuple_GET_ITEM(+args, 0), true};
        Object val{PyTuple_GET_ITEM(+args, 1), true};

        if (+key == SubClass<PyTypeObject>{&PyUnicode_Type}) {
            if (auto v = ref.request<Dictionary>()) {
                auto out = Object::from(PyDict_New());
                for (auto &x : *v) {
                    Object k = as_object(x.first);
                    Object v = python_cast(std::move(x.second), val, root);
                    if (!k || !v || PyDict_SetItem(out, k, v)) return {};
                }
                return out;
            }
        }

        if (auto v = ref.request<Vector<std::pair<Variable, Variable>>>()) {
            auto out = Object::from(PyDict_New());
            for (auto &x : *v) {
                Object k = python_cast(std::move(x.first), key, root);
                Object v = python_cast(std::move(x.second), val, root);
                if (!k || !v || PyDict_SetItem(out, k, v)) return {};
            }
            return out;
        }
    }
    return {};
}

// Convert Variable to a class which is a subclass of rebind.Variable
Object variable_cast(Variable &&v, Object const &t) {
    PyObject *x;
    if (t) x = +t;
    else if (!v.has_value()) return {Py_None, true};
    else if (auto it = python_types.find(v.type().info()); it != python_types.end()) x = +it->second;
    else x = type_object<Variable>();

    auto o = Object::from((x == type_object<Variable>()) ?
        PyObject_CallObject(x, nullptr) : PyObject_CallMethod(x, "__new__", "O", x));

    DUMP("making variable ", v.type());
    cast_object<Variable>(o) = std::move(v);
    DUMP("made variable ", v.type());
    return o;
}

/******************************************************************************/

Object bool_cast(Variable &&ref) {
    if (auto p = ref.request<bool>()) return as_object(*p);
    return {};
}

Object int_cast(Variable &&ref) {
    if (auto p = ref.request<Integer>()) return as_object(*p);
    DUMP("bad int");
    return {};
}

Object float_cast(Variable &&ref) {
    if (auto p = ref.request<double>()) return as_object(*p);
    if (auto p = ref.request<Integer>()) return as_object(*p);
    DUMP("bad float");
    return {};
}

Object str_cast(Variable &&ref) {
    DUMP("converting ", ref.type(), " to str");
    DUMP(ref.type(), bool(ref.action()));
    if (auto p = ref.request<std::string_view>()) return as_object(std::move(*p));
    if (auto p = ref.request<std::string>()) return as_object(std::move(*p));
    if (auto p = ref.request<std::wstring_view>())
        return {PyUnicode_FromWideChar(p->data(), static_cast<Py_ssize_t>(p->size())), false};
    if (auto p = ref.request<std::wstring>())
        return {PyUnicode_FromWideChar(p->data(), static_cast<Py_ssize_t>(p->size())), false};
    return {};
}

Object bytes_cast(Variable &&ref) {
    if (auto p = ref.request<BinaryView>()) return as_object(std::move(*p));
    if (auto p = ref.request<Binary>()) return as_object(std::move(*p));
    return {};
}

Object type_index_cast(Variable &&ref) {
    if (auto p = ref.request<TypeIndex>()) return as_object(std::move(*p));
    else return {};
}

Object function_cast(Variable &&ref) {
    if (auto p = ref.request<Function>()) return as_object(std::move(*p));
    else return {};
}

Object memoryview_cast(Variable &&ref, Object const &root) {
    if (auto p = ref.request<ArrayView>()) {
        auto x = type_object<ArrayBuffer>();
        auto obj = Object::from(PyObject_CallObject(x, nullptr));
        cast_object<ArrayBuffer>(obj) = {std::move(*p), root};
        return Object::from(PyMemoryView_FromObject(obj));
    } else return {};
}

Object getattr(PyObject *obj, char const *name) {
    if (PyObject_HasAttrString(obj, name))
        return {PyObject_GetAttrString(obj, name), false};
    return {};
}

// condition: PyType_CheckExact(type) is false
bool is_structured_type(PyObject *type, PyObject *origin) {
    if constexpr(PythonVersion >= Version(3, 7, 0)) {
        DUMP("is_structure_type 3.7A");
        // in this case, origin may or may not be a PyTypeObject *
        return origin == +getattr(type, "__origin__");
    } else {
        // case like typing.Union: type(typing.Union[int, float] == typing.Union)
        return (+type)->ob_type == reinterpret_cast<PyTypeObject *>(origin);
    }
}

bool is_structured_type(PyObject *type, PyTypeObject *origin) {
    if constexpr(PythonVersion >= Version(3, 7, 0)) {
        DUMP("is_structure_type 3.7B");
        return reinterpret_cast<PyObject *>(origin) == +getattr(type, "__origin__");
    } else {
        // case like typing.Tuple: issubclass(typing.Tuple[int, float], tuple)
        return is_subclass(reinterpret_cast<PyTypeObject *>(type), reinterpret_cast<PyTypeObject *>(origin));
    }
}

Object union_cast(Variable &&v, Object const &t, Object const &root) {
    if (auto args = type_args(t)) {
        auto const n = PyTuple_GET_SIZE(+args);
        for (Py_ssize_t i = 0; i != n; ++i) {
            Object o = python_cast(std::move(v), {PyTuple_GET_ITEM(+args, i), true}, root);
            if (o) return o;
            else PyErr_Clear();
        }
    }
    return type_error("cannot convert value to %R from type %S", +t, +type_index_cast(v.type()));
}


// Convert C++ Variable to a requested Python type
// First explicit types are checked:
// None, object, bool, int, float, str, bytes, TypeIndex, list, tuple, dict, Variable, Function, memoryview
// Then, the output_conversions map is queried for Python function callable with the Variable
Object try_python_cast(Variable &&v, Object const &t, Object const &root) {
    DUMP("try_python_cast ", v.type());
    if (auto it = type_translations.find(t); it != type_translations.end()) {
        DUMP("type_translation found");
        return try_python_cast(std::move(v), it->second, root);
    } else if (PyType_CheckExact(+t)) {
        auto type = reinterpret_cast<PyTypeObject *>(+t);
        DUMP("is Variable ", is_subclass(type, type_object<Variable>()));
        if (+type == Py_None->ob_type || +t == Py_None)       return {Py_None, true};                        // NoneType
        else if (type == &PyBool_Type)                        return bool_cast(std::move(v));                // bool
        else if (type == &PyLong_Type)                        return int_cast(std::move(v));                 // int
        else if (type == &PyFloat_Type)                       return float_cast(std::move(v));               // float
        else if (type == &PyUnicode_Type)                     return str_cast(std::move(v));                 // str
        else if (type == &PyBytes_Type)                       return bytes_cast(std::move(v));               // bytes
        else if (type == &PyBaseObject_Type)                  return as_deduced_object(std::move(v));        // object
        else if (is_subclass(type, type_object<Variable>()))  return variable_cast(std::move(v), t);         // Variable
        else if (type == type_object<TypeIndex>())            return type_index_cast(std::move(v));          // type(TypeIndex)
        else if (type == type_object<Function>())             return function_cast(std::move(v));            // Function
        else if (is_subclass(type, &PyFunction_Type))         return function_cast(std::move(v));            // Function
        else if (type == &PyMemoryView_Type)                  return memoryview_cast(std::move(v), root);    // memory_view
    } else {
        DUMP("Not type and not in translations");
        if (auto p = cast_if<TypeIndex>(t)) { // TypeIndex
            Dispatch msg;
            if (auto var = std::move(v).request_variable(msg, *p))
                return variable_cast(std::move(var));
            std::string c1 = v.type().name(), c2 = p->name();
            return type_error("could not convert object of type %s to type %s", c1.data(), c2.data());
        }
        else if (is_structured_type(t, UnionType))     return union_cast(std::move(v), t, root);
        else if (is_structured_type(t, &PyList_Type))  return list_cast(std::move(v), t, root);       // List[T] for some T (compound type)
        else if (is_structured_type(t, &PyTuple_Type)) return tuple_cast(std::move(v), t, root);      // Tuple[Ts...] for some Ts... (compound type)
        else if (is_structured_type(t, &PyDict_Type))  return dict_cast(std::move(v), t, root);       // Dict[K, V] for some K, V (compound type)
        DUMP("Not one of the structure types");
    }

    DUMP("custom convert ", output_conversions.size());
    if (auto p = output_conversions.find(t); p != output_conversions.end()) {
        DUMP(" conversion ");
        Object o = variable_cast(std::move(v));
        if (!o) return type_error("could not cast Variable to Python object");
        DUMP("calling function");
        auto &obj = static_cast<Var &>(cast_object<Variable>(o)).ward;
        if (!obj) obj = root;
        return Object::from(PyObject_CallFunctionObjArgs(+p->second, +o, nullptr));
    }

    return nullptr;
}

Object python_cast(Variable &&v, Object const &t, Object const &root) {
    Object out = try_python_cast(std::move(v), t, root);
    if (!out) return type_error("cannot convert value to type %R from type %S", +t, +type_index_cast(v.type()));
    return out;
}

/******************************************************************************/

}
