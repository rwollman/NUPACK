import inspect, importlib, enum, typing, sys

################################################################################

class Scalar(enum.IntEnum):
    Bool = 0
    Char = 1
    SignedChar = 2
    UnsignedChar = 3
    Unsigned = 4
    Signed = 5
    Float = 6
    Pointer = 7

################################################################################

scalar_types = {
    'float32': (Scalar.Float, 32),
    'float64': (Scalar.Float, 64),
    'uint16':  (Scalar.Unsigned, 16),
    'uint32':  (Scalar.Unsigned, 32),
    'uint64':  (Scalar.Unsigned, 64),
    'int16':   (Scalar.Signed, 16),
    'int32':   (Scalar.Signed, 32),
    'int64':   (Scalar.Signed, 64)
}

def find_scalars(info):
    find = lambda x, n: next(t for c, t, s in info if c == x and s == n)
    return {k: find(x, n) for k, (x, n) in scalar_types.items()}

################################################################################

translations = {
    '{}': '__str__',
    '[]': '__getitem__',
    '()': '__call__',
    '~': '__invert__',
    '^': '__xor__',
    '+': '__add__',
    '-': '__sub__',
    '/': '__div__',
    '*': '__mul__',
    '==': '__eq__',
    '!=': '__ne__',
    '<': '__lt__',
    '>': '__gt__',
    '<=': '__le__',
    '>=': '__ge__',
    'bool': '__bool__',
}

################################################################################

def unwrap(old):
    return getattr(old, '_origin_', old)

 # The following functions are just used for their signatures

def opaque_signature(*args):
    '''Call a function from a rebind.Document'''

def default_str(self) -> str:
    '''Convert self to a str via C++'''

def default_bool(self) -> bool:
    '''Convert self to a boolean via C++'''

def default_logical(self, other, _fun_=None):
    '''Run a logical operation via C++'''
    if type(self) != type(other):
        return NotImplemented
    return _fun_(self, other).cast(bool)

def default_int(self) -> int:
    '''Run an integer operation via C++'''

default_methods = {'__%s__' % k: f for k, f in (
    ('str', default_str),
    ('repr', default_str),
    ('bool', default_bool),
    ('contains', default_logical),
    *[(k, default_logical) for k in 'eq ne lt gt le ge'.split()],
    *[(k, default_int) for k in 'int len index hash'.split()]
)}

################################################################################

def signatures(function):
    '''Get signatures out of a Function or a wrapped Function'''
    try:
        function = function.__kwdefaults__['_orig']
    except AttributeError:
        pass
    return function.signatures()

################################################################################

def finalize(func, log):
    '''Finalize C++ held Python objects from rebind'''
    log.info('cleaning up Python C++ resources')
    func()

################################################################################

def import_namespace(name):
    try:
        return importlib.import_module(name)
    except ImportError:
        x, y = name.rsplit('.', maxsplit=1)
        return getattr(import_namespace(x), y)

################################################################################

def split_module(pkg, name):
    '''Return the module containing pkg.name, and the last string in name'''
    x, y = '{}.{}'.format(pkg, name).rsplit('.', maxsplit=1)
    return import_namespace(x), y

################################################################################

def discard_parameter(sig, key):
    return inspect.Signature([v for v in sig.parameters.values() if v.name != key],
            return_annotation=sig.return_annotation)

################################################################################

def update_module(dest, source, pattern):
    if isinstance(source, str):
        source = importlib.import_module(source)
    if pattern == '*':
        pattern = getattr(source, '__all__', pattern)
    if pattern == '*':
        for k in dir(source):
            if not k.startswith('_'):
                setattr(dest, k, getattr(source, k))
    else:
        for k in (pattern.split() if isinstance(pattern, str) else pattern):
            setattr(dest, k, getattr(source, k))
