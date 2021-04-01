namespace rebind {

// move_from is called 1) during init, V.move_from(V), to transfer the object (here just use Var move constructor)
//                     2) during assignment, R.move_from(L), to transfer the object (here cast V to new object of same type, swap)
//                     2) during assignment, R.move_from(V), to transfer the object (here cast V to new object of same type, swap)
PyObject * var_copy_assign(PyObject *self, PyObject *value) noexcept {
    return raw_object([=] {
        DUMP("- copying variable");
        cast_object<Var>(self).assign(variable_reference_from_object({value, true}));
        return Object(self, true);
    });
}

PyObject * var_move_assign(PyObject *self, PyObject *value) noexcept {
    return raw_object([=] {
        DUMP("- moving variable");
        auto &s = cast_object<Var>(self);
        Variable v = variable_reference_from_object({value, true});
        v.move_if_lvalue();
        s.assign(std::move(v));
        // if (auto p = cast_if<Variable>(value)) p->reset();
        return Object(self, true);
    });
}

/******************************************************************************/

int var_bool(PyObject *self) noexcept {
    if (auto v = cast_if<Variable>(self)) return v->has_value();
    else return PyObject_IsTrue(self);
}

PyObject * var_has_value(PyObject *self, PyObject *) noexcept {
    return PyLong_FromLong(var_bool(self));
}

/******************************************************************************/

PyObject * var_cast(PyObject *self, PyObject *type) noexcept {
    return raw_object([=] {
        return python_cast(std::move(cast_object<Variable>(self)), Object(type, true), Object(self, true));
    });
}

PyObject * var_type(PyObject *self, PyObject *) noexcept {
    return raw_object([=] {
        auto o = Object::from(PyObject_CallObject(type_object<TypeIndex>(), nullptr));
        cast_object<TypeIndex>(o) = cast_object<Variable>(self).type();
        return o;
    });
}

PyObject * var_qualifier(PyObject *self, PyObject *) noexcept {
    return raw_object([=] {
        return as_object(static_cast<Integer>(cast_object<Variable>(self).qualifier()));
    });
}

PyObject * var_is_stack_type(PyObject *self, PyObject *) noexcept {
    return raw_object([=] {
        return as_object(cast_object<Variable>(self).is_stack_type());
    });
}

PyObject * var_address(PyObject *self, PyObject *) noexcept {
    return raw_object([=] {
        return as_object(Integer(reinterpret_cast<std::uintptr_t>(cast_object<Variable>(self).data())));
    });
}

PyObject * var_ward(PyObject *self, PyObject *) noexcept {
    return raw_object([=] {
        Object out = cast_object<Var>(self).ward;
        return out ? out : Object(Py_None, true);
    });
}

PyObject * var_set_ward(PyObject *self, PyObject *arg) noexcept {
    return raw_object([=]() -> Object {
        Object root{arg, true};
        while (true) { // recurse upwards to find the governing lifetime
            auto p = cast_if<Var>(root);
            if (!p || !p->ward) break;
            root = p->ward;
        }
        cast_object<Var>(self).ward = std::move(root);
        return {self, true};
    });
}

PyObject * var_from(PyObject *cls, PyObject *obj) noexcept {
    return raw_object([=]() -> Object {
        if (PyObject_TypeCheck(obj, reinterpret_cast<PyTypeObject *>(cls))) {
            // if already correct type
            return {obj, true};
        }
        if (auto p = cast_if<Variable>(obj)) {
            // if a Variable try .cast
            return python_cast(static_cast<Variable const &>(*p).reference(), Object(cls, true), Object(obj, true));
        }
        // Try cls.__init__(obj)
        return Object::from(PyObject_CallFunctionObjArgs(cls, obj, nullptr));
    });
}

/******************************************************************************/

PyNumberMethods VarNumberMethods = {
    .nb_bool = static_cast<inquiry>(var_bool),
};

PyMethodDef VarMethods[] = {
    {"copy_from",     static_cast<PyCFunction>(var_copy_assign),   METH_O,       "assign from other using C++ copy assignment"},
    {"move_from",     static_cast<PyCFunction>(var_copy_assign),   METH_O,       "assign from other using C++ move assignment"},
    {"address",       static_cast<PyCFunction>(var_address),       METH_NOARGS,  "get C++ pointer address"},
    {"_ward",         static_cast<PyCFunction>(var_ward),          METH_NOARGS,  "get ward object"},
    {"_set_ward",     static_cast<PyCFunction>(var_set_ward),      METH_O,       "set ward object and return self"},
    {"qualifier",     static_cast<PyCFunction>(var_qualifier),     METH_NOARGS,  "return qualifier of self"},
    {"is_stack_type", static_cast<PyCFunction>(var_is_stack_type), METH_NOARGS,  "return if object is held in stack storage"},
    {"type",          static_cast<PyCFunction>(var_type),          METH_NOARGS,  "return TypeIndex of the held C++ object"},
    {"has_value",     static_cast<PyCFunction>(var_has_value),     METH_NOARGS,  "return if a C++ object is being held"},
    {"cast",          static_cast<PyCFunction>(var_cast),          METH_O,       "cast to a given Python type"},
    {"from_object",   static_cast<PyCFunction>(var_from),          METH_CLASS | METH_O, "cast an object to a given Python type"},
    {nullptr, nullptr, 0, nullptr}
};

template <>
PyTypeObject Holder<Var>::type = []{
    auto o = type_definition<Var>("rebind.Variable", "C++ class object");
    o.tp_as_number = &VarNumberMethods;
    o.tp_methods = VarMethods;
    // no init (just use default constructor)
    // tp_traverse, tp_clear
    // PyMemberDef, tp_members
    return o;
}();

/******************************************************************************/

}
