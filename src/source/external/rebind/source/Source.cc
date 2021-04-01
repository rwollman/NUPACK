#include <rebind/Document.h>

/******************************************************************************/

#if __has_include(<cxxabi.h>)
#   include <cxxabi.h>
    namespace rebind::runtime {
        using namespace __cxxabiv1;

        std::string demangle(char const *s) {
            int status = 0;
            char * buff = __cxa_demangle(s, nullptr, nullptr, &status);
            if (!buff) return s;
            std::string out = buff;
            std::free(buff);
            return out;
        }

        char const *unknown_exception_description() noexcept {
            return abi::__cxa_current_exception_type()->name();
        }
    }
#else
    namespace rebind::runtime {
        std::string demangle(char const *s) {return s;}

        char const *unknown_exception_description() noexcept {return "C++: unknown exception";}
    }
#endif

/******************************************************************************/

namespace rebind {

Demangler demangler{&runtime::demangle};

void set_demangler(Demangler fun) noexcept {demangler = std::move(fun);}

std::string demangle(char const *s) {
    if (demangler) return demangler(s);
    else return s;
}

/******************************************************************************/

bool Debug = false;

void set_debug(bool debug) noexcept {Debug = debug;}
bool debug() noexcept {return Debug;}

/******************************************************************************/

Document & document() noexcept {
    static Document static_document;
    return static_document;
}

void render_default(Document &, std::type_info const &t) {
    if (debug()) std::cout << "Not rendering type " << t.name() << std::endl;
}

/******************************************************************************/

void lvalue_fails(Variable const &v, Dispatch &msg, TypeIndex t) {
    char const *s = "could not convert to lvalue reference";
    if (v.type() == t) {
        if (v.qualifier() == Rvalue) s = "could not convert rvalue to lvalue reference";
        if (v.qualifier() == Const) s = "could not convert const value to lvalue reference";
        if (v.qualifier() == Value) s = "could not convert value to lvalue reference";
    }
    msg.error(s, t);
}
/******************************************************************************/

void rvalue_fails(Variable const &v, Dispatch &msg, TypeIndex t) {
    char const *s = "could not convert to rvalue reference";
    if (v.type() == t) {
        if (v.qualifier() == Lvalue) s = "could not convert lvalue to rvalue reference";
        if (v.qualifier() == Const) s = "could not convert const value to rvalue reference";
    }
    msg.error(s, t);
}
/******************************************************************************/

void Variable::assign(Variable v) {
    // if *this is a value, make a copy or move the right hand side in
    if (qualifier() == Value) {
        if (v.qualifier() == Value) {
            *this = std::move(v);
        } else {
            DUMP("assign1");
            if (auto p = handle())
                act(ActionType::destroy, pointer(), nullptr);
            // Copy data but set the qualifier to Value
            static_cast<VariableData &>(*this) = v;
            idx.set_qualifier(Value);
            // DUMP(name(), qualifier(), name(), stack, v.stack);
            // DUMP(&v, v.pointer());
            // e.g. value = lvalue which means
            // Move variable if it held RValue
            act((v.qualifier() == Rvalue) ? ActionType::move : ActionType::copy, v.pointer(), this);
            // DUMP(stack, v.stack);
        }
    } else if (qualifier() == Const) {
        DUMP("assign bad");
        throw std::invalid_argument("Cannot assign to const Variable");
    } else { // qual == Lvalue or Rvalue
        DUMP("assigning reference", type(), &buff, pointer(), v.type());
        // qual, type, etc are unchanged
        act(ActionType::assign, pointer(), &v);
        if (v.has_value())
            throw std::invalid_argument("Could not coerce Variable to matching type");
    }
}

/******************************************************************************/

Variable Variable::request_var(Dispatch &msg, TypeIndex const &t, Qualifier q) const {
    DUMP((act != nullptr), " asking for ", t, "from", q, type());
    Variable v;
    if (!has_value()) {
        // Nothing to do; request always fails
    } else if (idx.matches(t)) { // Exact type match
        // auto info = reinterpret_cast<std::type_info const * const &>(t);
        if (t.qualifier() == Value) { // Make a copy or move
            v.idx = t;
            v.act = act;
            v.stack = stack;
            act((q == Rvalue) ? ActionType::move : ActionType::copy, pointer(), &v);
        } else if (t.qualifier() == Const || t.qualifier() == q) { // Bind a reference
            DUMP("yope", t, type(), q);
            reinterpret_cast<void *&>(v.buff) = pointer();
            v.idx = t;
            v.act = act;
            v.stack = stack;
        } else {
            DUMP("nope");
            msg.error("Source and target qualifiers are not compatible");
        }
    } else {
        ::new(static_cast<void *>(&v.buff)) RequestData{t, &msg, q};
        act(ActionType::response, pointer(), &v);
        // DUMP(v.has_value(), v.name(), q, v.qualifier());

        if (!v.has_value()) {
            DUMP("response returned no value");
            msg.error("Did not respond with anything");
        } else if (v.type() != t)  {
            DUMP("response gave wrong type", v.type(), t);
            msg.error("Did not respond with correct type");
            v.reset();
        }
    }
    DUMP(v.type(), t);
    return v;
}

/******************************************************************************/

void set_source(Dispatch &msg, std::type_info const &t, Variable &&v) {
    if (auto p = std::move(v).target<std::string &&>()) {
        msg.source = std::move(*p);
    } else if (auto p = v.target<std::string_view const &>()) {
        msg.source = *p;
    } else if (auto p = v.target<TypeIndex const &>()) {
        msg.source = p->name();
    } else {
        msg.source = t.name();
    }
}

/******************************************************************************/

TypeData & Document::type(TypeIndex t, std::string s, Variable data) {
    auto it = contents.try_emplace(std::move(s), TypeData()).first;
    if (auto p = it->second.target<TypeData &>()) {
        p->data[t] = std::move(data);
        types[t] = &(*it);
        return *p;
    }
    DUMP(type_index<decltype(it->second)>(), " ", t, " ", s, " ", data.type(), " ", it->second.type());
    s.insert(0, "tried to declare both a non-type and a type for the same key ");
    throw std::runtime_error(std::move(s));
}

Function & Document::find_method(TypeIndex t, std::string name) {
    if (auto it = types.find(t); it != types.end()) {
        if (auto p = it->second->second.target<TypeData &>())
            return p->methods.emplace(std::move(name), Function()).first->second;
        name.insert(0, "tried to declare a method ");
        name += "for a non-type key ";
        name += it->second->first;
        throw std::runtime_error(std::move(name));
    }
    name.insert(0, "tried to declare a method ");
    name += "for the undeclared type ";
    name += t.name();
    throw std::runtime_error(std::move(name));
}

Function & Document::find_function(std::string s) {
    auto it = contents.emplace(std::move(s), Type<Function>()).first;
    if (auto f = it->second.target<Function &>()) return *f;
    throw std::runtime_error("tried to declare both a non-function and a function for the same key " + it->first);
}

/******************************************************************************/

}
