# Python API

The `rebind` Python API is generally designed around the following workflow:

1. A C++ API is defined for your classes and functions. It's exported into a module object like `my_cpp_module.so` into your Python package directory next to `__init__.py`.
2. A placeholder Python API is defined for your whole module, where you can add customizations and documentation in the `.py` files.
3. When importing your module, the `rebind` Python module is used to fill out all of the placeholders with their C++ implementations.

Here's an example of a pretty minimal `__init__.py` that implements Step 3 above.

```python
# import whatever modules you want as normal in your __init__.py here

try:
    from .my_cpp_module import document
except ImportError:
    raise ImportError('C++ module object "my_package.my_cpp_module" cannot be imported')

from rebind import render_module
rendered_document, config = render_module(__name__, document)

# if you need to do some additional steps after export, you can do them here (rare)
```

In this wrapper code, there are 3 extra artifacts of the export process. They can all be discarded if you want.

1. `document` is the raw document listing all of the functions and types that are exported from your C++ module.
2. `rendered_document` is the modified document which includes maps from C++ entities to their Python equivalents.
3. `config` is an instance of `rebind.Config` which can be used to mess with global `rebind` functionality for your module.

## Python Config

The `rebind.Config` class is pretty simple and offers the following methods for you to use:

1. `set_input_conversion` sets an automatic conversion from one Python type to another. For instance, if you use `numpy` but need all of your input arrays to be column-major once they hit C++, you could add this conversion:
```python
config.set_input_conversion(numpy.ndarray, numpy.asfortranarray)
```
1. `set_output_conversion` adds a custom method to accomplish conversion of a `rebind.Variable` into a desired Python class. For instance, `rebind` doesn't provide a direct `numpy` interface, but you can still return a contiguous array by using a `Variable -> memoryview -> numpy.ndarray` conversion sequence:
```python
config.set_output_conversion(numpy.ndarray, lambda variable: numpy.asarray(variable.cast(memoryview)))
```
3. `debug` is an instance property with get/set methods to turn on `rebind` printing debug messages to `stdout`:
```python
config.debug = True
```

## Wrapping a C++ function

### Exporting a function without writing Python code

```c++
doc.function("mymodule.add_float_to_int", [](int x, double y) {
    return x + y;
});
```

If you write no corresponding Python code in `mymodule`, `rebind` will define a function with no docstring, argument names, or default arguments. (This function will accept positional arguments only.) In addition, this function will return the native C++ type, which in this case is an instance of `rebind.Float`:

```python
mymodule.add_float_to_int(1, 2.5) # -> Float(3.5)
```

`rebind` defines its own builtin types to more closely mirror their C++ implementation. For instance, all of these types are mutable:

```python
z = mymodule.add_float_to_int(1, 2.5)
z2 = z
z2 += 1
print(z) # -> Float(4.5)
```

This provides the most flexible raw approach. However, usually you will want to wrap the function using Python code as is described below.

### Creating an annotated function

Let's say we want to customize the behavior of `add_float_to_int` in 4 ways:

1. Add a docstring to the function.
2. Add keyword names to the function.
3. Add argument defaults to the function.
4. Automatically cast the result to a Python `float`.

We can do all of this by writing the simple wrapper code below in `mymodule.py`:

```python
def add_float_to_int(x, y=0.0) -> float:
    '''Adds float to an int'''
```

All this uses the Python `rebind` API to create a nice syntactic sugar for what we want. Now we can use the code expectedly as follows:

```python
z = mymodule.add_float_to_int(x=5)
print(z, type(z)) # -> 5, float
```

Note that if we left off the return annotation, we could still return a `rebind.Float` if we wanted to.

### Overriding a function with completely custom behavior

Sometimes, we want to invoke more complicated wrapping functionality than described above. In these cases, you can use the `_fun_` based API. To use this API, define your function with a defaulted keyword argument `_fun_=None`. A function consumer will never use this keyword argument directly. Instead, `rebind` wil set `_fun_` to the original raw exported function for your wrapper code in the function body. For example, here's a contrived example of adding some more complicated functionality to `add_float_to_int`.

```python
def add_float_to_int(x, y=None, _fun_=None):
    '''Adds float to an int'''
    if y is None:
        y = 5.5
    elif not isinstance(y, float):
        raise TypeError('Expected y to be a float but got {}'.format(type(y)))
    return _fun_(x, y).cast(float) # Note that return type casting must be done in the function body with this API
```

### Additional details

For the most part, you are intended to pass only positional arguments to `rebind.Function`. However, the call operator also accepts the keywords `gil`, `return_type`, and `signature`. These are for more custom behavior as described below.

#### Manipulating the GIL

If specified, `gil` is expected to be a `bool` of whether the Python global interpreter lock should be held (`True`) or released (`False`). `gil` defaults to `True`. However, for long-running `C++` code, you can turn off the `gil` to allow multiple simultaneous threads to execute.

Note that if you specify `gil=False` but call a Python callback from your C++ code, `rebind` will automatically re-acquire the GIL during the scope of the callback. That means you shouldn't have to worry about segfaulting in any case. The general reason to leave `gil=True` is to avoid overhead for simple functions or functions that are not expected to execute concurrently.

#### Manually choosing an overload

The `rebind` approach to overloading is generally to try one overload after another until one works. However, you might want to circumvent this process for performance or another reason. To help `rebind` to choose the correct overload, you can specify *either* `return_type` or `signature`.

Specifying `return_type` is a relatively simple way to dispatch on the returned type of the overload. It's especially useful for use inside `__init__` methods to choose a specific type to instantiate from Python.


## Wrapping a C++ class

Whenever you export a C++ class, you'll generally want to define a corresponding Python one. Suppose we have this class:

```c++
struct Adder {
    bool add = false;
    double operator()(int x, double y=5) const {
        return add ? x + y : x;
    }
};
```

And we export it as follows:

```
Type<Adder> t;
doc.type("mymodule.Adder", t);
doc.method(t, "new", rebind::construct<bool>(t));
doc.method(t, "()", &Adder::operator());
doc.member(t, &Adder::add);
```

A minimal wrapper will just define the class:

```python
class Adder:
    pass
```

Actually, things can get a bit confusing if you have messed up your export name. For example, if you export `Adder` to the wrong module, you might just end up with a blank `Adder` class like written here. Therefore it's advised to prevent any unintened instantiation via a disabling mechanism like this:

```python
class Adder:
    __new__ = NotImplemented
```

This is convenient because you will generally never override the `__new__` method of a `rebind` class. It's not necessary, but it could save you some headaches when debugging your code.

### Class members

A more complete wrapper might be written like below

```python
class Adder:
    ```An example class```
    add: bool

    __new__ = NotImplemented

    def __init__(self, add=False):
        '''Initialize based on whether adding should occur'''

    def __call__(self, x, y=5) -> float:
        '''Either add x and y or just return x'''
```

These method wrappers work pretty much the same as function wrappers in `rebind`, and the docstring part is pretty obvious. The main wrinkle is on exporting a class instance member via `add: bool`. You don't have to do this, in which case the member will still be exported, but providing a type annotation will make `Adder().add` look like a Python `bool`. In addition, annotating the class members is a good practice for documenting the class in any case. If you really want to document the member without assigning a type, you can use `add: None`.

## Wrapping a global or static variable
