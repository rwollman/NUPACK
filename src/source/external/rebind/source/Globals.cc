#include <rebind-python/API.h>

namespace rebind {

Object UnionType, TypeError;

std::unordered_map<Object, Object> type_translations{}, output_conversions{}, input_conversions{};

std::unordered_map<std::type_index, Object> python_types{};


void initialize_global_objects() {
    TypeError = {PyExc_TypeError, true};

    auto t = Object::from(PyImport_ImportModule("typing"));
    UnionType = Object::from(PyObject_GetAttrString(t, "Union"));
    // (+u)->ob_type
}

void clear_global_objects() {
    input_conversions.clear();
    output_conversions.clear();
    type_translations.clear();
    python_types.clear();
    UnionType = nullptr;
    TypeError = nullptr;
}

std::unordered_map<TypeIndex, std::string> type_names = {
    {typeid(void),             "void"},
    {typeid(void *),           "pointer"},
    {typeid(PyObject),         "PyObject"},
    {typeid(PyObject *),       "PyObject *"},
    {typeid(bool),             "bool"},
    {typeid(Real),             "float64"},
    {typeid(std::string_view), "str"},
    {typeid(std::string),      "str"},
    {typeid(TypeIndex),        "TypeIndex"},
    {typeid(Binary),           "Binary"},
    {typeid(BinaryView),       "BinaryView"},
    {typeid(BinaryData),       "BinaryData"},
    {typeid(ArrayView),        "ArrayView"},
    {typeid(Function),         "Function"},
    {typeid(Variable),         "Variable"},
    {typeid(Sequence),         "Sequence"},
    {typeid(char),             "char"},
    {typeid(unsigned char),    "unsigned_char"},
    {typeid(signed char),      "signed_char"},
    {typeid(char16_t),         "char16_t"},
    {typeid(char32_t),         "char32_t"},
    {typeid(int),              "int32"},
    {typeid(float),            "float32"},
    {typeid(long double),      "long_double"},
    {typeid(std::uint8_t),     "uint8"},
    {typeid(std::uint16_t),    "uint16"},
    {typeid(std::uint32_t),    "uint32"},
    {typeid(std::uint64_t),    "uint64"},
    {typeid(std::int8_t),      "int8"},
    {typeid(std::int16_t),     "int16"},
    {typeid(std::int32_t),     "int32"},
    {typeid(std::int64_t),     "int64"}
};


Zip<std::string_view, std::type_info const *> Buffer::formats = {
    {"d", &typeid(double)},
    {"f", &typeid(float)},
    {"c", &typeid(char)},
    {"b", &typeid(signed char)},
    {"B", &typeid(unsigned char)},
    {"?", &typeid(bool)},
    {"h", &typeid(short)},
    {"H", &typeid(unsigned short)},
    {"i", &typeid(int)},
    {"I", &typeid(unsigned int)},
    {"l", &typeid(long)},
    {"L", &typeid(unsigned long)},
    {"q", &typeid(long long)},
    {"Q", &typeid(unsigned long long)},
    {"n", &typeid(ssize_t)},
    {"s", &typeid(char[])},
    {"p", &typeid(char[])},
    {"N", &typeid(size_t)},
    {"P", &typeid(void *)}
};

#define REBIND_TMP(C, T) {Scalar::C, typeid(T), sizeof(T) * CHAR_BIT}

Zip<Scalar, TypeIndex, unsigned> scalars = {
    REBIND_TMP(Bool,         bool),
    REBIND_TMP(Char,         char),
    REBIND_TMP(SignedChar,   signed char),

    REBIND_TMP(UnsignedChar, unsigned char),
    REBIND_TMP(UnsignedChar, char16_t),
    REBIND_TMP(UnsignedChar, char32_t),

    REBIND_TMP(Unsigned,     std::uint8_t),
    REBIND_TMP(Unsigned,     std::uint16_t),
    REBIND_TMP(Unsigned,     std::uint32_t),
    REBIND_TMP(Unsigned,     std::uint64_t),

    REBIND_TMP(Unsigned,     unsigned short),
    REBIND_TMP(Unsigned,     unsigned int),
    REBIND_TMP(Unsigned,     unsigned long),
    REBIND_TMP(Unsigned,     unsigned long long),

    REBIND_TMP(Signed,       std::int8_t),
    REBIND_TMP(Signed,       std::int16_t),
    REBIND_TMP(Signed,       std::int32_t),
    REBIND_TMP(Signed,       std::int64_t),

    REBIND_TMP(Signed,       short),
    REBIND_TMP(Signed,       int),
    REBIND_TMP(Signed,       long),
    REBIND_TMP(Signed,       long long),

    REBIND_TMP(Float,        float),
    REBIND_TMP(Float,        double),
    REBIND_TMP(Float,        long double),
    REBIND_TMP(Pointer,      void *),
};
#undef REBIND_TMP

}
