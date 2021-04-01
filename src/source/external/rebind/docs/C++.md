
## User conversions

### Request

`Request` converts a `Variable` to type `T`. It returns a null object if this is impossible. To specify more information on an error, call `msg.error(...)`.

In this case `T` may be a reference type. In the below, read `return_type` as:
- `T *` for `Request<T &>`
- `T *` for `Request<T &&>`
- `T const *` for `Request<T const &>`
- `std::optional<T>` for `Request<T>`


```c++
template <class T, class SFINAE=void>
struct Request {
    /* Mandatory */ return_type operator()(Variable const &r, Dispatch &msg);
    /* Optional  */ return_type operator()(Variable &&r, Dispatch &msg);
};
```

### Response

`Response<T>` converts type `T` into a `Variable` holding the requested type `idx`. If the conversion is impossible, it does nothing. In this case, `T` is always a non-reference type.

```c++
template <class T, class SFINAE=void>
struct Response {
    /* Mandatory */ void operator()(Variable &out, TypeIndex idx, T const &t);
    /* Optional  */ void operator()(Variable &out, TypeIndex idx, T &&t);
};
```

### Renderer

`Renderer` modifies a `Document`.

```c++
template <class T, class SFINAE=void>
struct Renderer {
    void operator()(Document &doc);
};
```

## Variable

### Constructors

```c++
constexpr Variable() noexcept = default;
```

Construct from reference type
```c++
// Reference type
template <class T, std::enable_if_t<!(std::is_same_v<std::decay_t<T>, T>), int> = 0>
Variable(Type<T>, typename SameType<T>::type t) noexcept;
```

Construct from non-reference type
```c++
template <class T, class ...Ts, std::enable_if_t<(std::is_same_v<std::decay_t<T>, T>), int> = 0>
Variable(Type<T>, Ts &&...ts);

template <class T, std::enable_if_t<!std::is_base_of_v<VariableData, unqualified<T>>, int> = 0>
Variable(T &&t);
```

Copy, move, destruct
```c++
Variable(Variable &&v) noexcept;
Variable(Variable const &v);
Variable & operator=(Variable &&v) noexcept;
Variable & operator=(Variable const &v);
~Variable();
```

### Mutations

```c++
void reset();
void assign(Variable v);
```

### Descriptors

```c++
void const * data() const;
char const * name() const;
std::type_index type() const;
std::type_info const & info() const;
Qualifier qualifier() const;
ActionFunction action() const;
bool is_stack_type() const;

constexpr bool has_value() const;
explicit constexpr operator bool() const;
```

### Make a reference
```c++
Variable reference() &;
Variable reference() const &;
Variable reference() &&;
```

### Conversions

Request any type T by non-custom conversions
```c++
Variable request_variable(Dispatch &msg, std::type_index const, Qualifier q=Value) const;
```

Request reference T by custom conversions
```c++
template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
std::remove_reference_t<T> *request(Dispatch &msg={}, Type<T> t={}) const;

template <class T, std::enable_if_t<!std::is_reference_v<T>, int> = 0>
std::optional<T> request(Dispatch &msg={}, Type<T> t={}) const;
```

Request non-reference T by custom conversions
```c++
template <class T>
T downcast(Dispatch &msg={}, Type<T> t={}) const;
```

Return pointer to target if it is trivially convertible to requested type
```c++
template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
std::remove_reference_t<T> *target(Type<T> t={}) const &;

```
Return pointer to target if it is trivially convertible to requested type
```c++
template <class T, std::enable_if_t<std::is_reference_v<T>, int> = 0>
std::remove_reference_t<T> *target(Type<T> t={}) &;
```

## Document

In C++, without necessary python inclusion, we define a document which more or less represents a module.

It is a map from names to free functions and classes.

- Document:
    - map of `String` to (`Function` or `Class`)
- Function: a simple immutable callable object
    - `Value(ArgPack)`
    - `Vector<String>` keywords (also types?)
- Class: an object containing named members and methods
    - map of `String` to (`Member` or `Method`)
- Member: provides access to a member given an instance of a C++ class
    - summarizable as `Value(Any &)`
- Value: a C++ type which is a variant of
    - `Any`
    - `Real`
    - `Integer`
    - ... other built-ins
- Method: same as function but passes `this` by reference as first argument
    - `Value(Any &, ArgPack)`
    - `Vector<String> keywords`

### Members

```python
variable.member # property that returns a reference to the member with parent as a ward
variable.member = other_variable # setter that calls copy assign from other to member
variable.member.move_from(other_variable) # setter that calls move assign from other to member
variable.move_from(other_variable) # if variable is V, move_from,
```

## List of good pybind11 features

- possibly pypy
- static fields, methods (maybe)
- dynamic attributes

## Summary

These two sides are kept pretty modular; in particular, your C++ testing code doesn't need to know anything about (or include) any Python headers. In terms of API, `rebind` draws on ideas from `Catch` and `doctest` but tries to be less macro-based.

There are a lot of nice features in `rebind` including the abilities to:
- write test event handlers and CLIs in Python (easier than in C++)
- run tests in a threadsafe and parallel manner
- parametrize your tests either in C++, Python, or from the command line
- call tests from other tests within C++ or Python
- replace the given Python handlers with your own without any recompilation of your C++ code

The primary hurdles to using `rebind` are that it requires C++17 (most importantly for `std::variant`) and Python 2.7+/3.3+ for the Python reporters. The build can also be bit more tricky, though it's not as complicated as you might expect.

`rebind` is quite usable but also at a pretty early stage of development. Please try it out and give some feedback!

### Canonical codes
- `new`: constructor
- `()`: class call operator

### Reference semantics
As the caller we can call our function with
- `T &`: we expect that our object can be mutated
- `T const &`: we do not expect that our object can be mutated
- `T &&`: we give up our object and don't care about it

As a function we may wish to be called with
- `T`: we don't want to modify the caller, but we do want to modify the input
- `T &`: we do want to modify both
- `T const &`: we don't want to modify either
- `T &&`: we take it

We can convert as follows:
- `T`: we can always convert to this
- `T &`: take `T &` or else make a copy and pass that instead
- `T const &`: we don't want to modify either
- `T &&`: take `T &&` or else make a copy and pass that instead
