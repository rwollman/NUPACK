import inspect, importlib, functools, logging, typing, atexit, collections
from . import Config, common, ConversionError

log = logging.getLogger(__name__)
# logging.basicConfig(level=logging.INFO)

################################################################################

def render_module(pkg: str, doc: dict, set_type_names=False, monkey=[]):
    clear = doc['clear_global_objects']
    try:
        return _render_module(pkg, doc, set_type_names, monkey=monkey)
    except BaseException:
        clear()
        raise
    finally:
        atexit.register(common.finalize, clear, log)

################################################################################

def _render_module(pkg, doc, set_type_names, monkey=[]):
    log.info('rendering document into module %s', repr(pkg))
    config, out = Config(doc), doc.copy()
    config.set_type_error(ConversionError)

    classes, modules, translate = set(), set(), {}

    # render classes and methods
    out['types'] = {k: v for k, v in doc['contents'] if isinstance(v, tuple)}
    for k, (meth, data) in out['types'].items():
        mod, cls = render_type(translate, pkg, (doc['Variable'],), k, dict(meth))
        modules.add(mod)
        classes.add(cls)
        cls._metadata_ = {k: v or None for k, v in data}
        for k, v in data:
            config.set_type(k, cls)

    # put type names if desired
    if set_type_names:
        names = tuple((o[0], k) for k, v in out['types'].items() for o in v[1])
        config.set_type_names(names)

    # render global objects (including free functions)
    obj = lambda old, new: translate.__setitem__(old, new) or new
    out['objects'] = {k: obj(*render_object(pkg, k, v)) for k, v in doc['contents'] if not isinstance(v, tuple)}

    out['scalars'] = common.find_scalars(doc['scalars'])

    # Monkey-patch modules based on redefined types (takes care of simple cases at least)
    mod = importlib.import_module(pkg)

    for mod in tuple(modules):
        parts = mod.__name__.split('.')
        for i in range(1, len(parts)):
            try:
                modules.add(importlib.import_module('.'.join(parts[:i])))
            except ImportError:
                pass

    translate.pop(None, None)

    for mod in modules.union(classes).union(monkey):
        log.info('rendering monkey-patching namespace {}'.format(mod))
        for k in dir(mod):
            try:
                setattr(mod, k, translate[getattr(mod, k)])
                log.info('monkey-patching namespace {} attribute {}'.format(mod, k))
            except (TypeError, KeyError):
                pass

    for k, v in translate.items():
        if isinstance(k, type) and isinstance(v, type):
            config.set_translation(k, v)

    log.info('finished rendering document into module %s', repr(pkg))
    return out, config

################################################################################

def render_init(init):
    if init is None:
        def __init__(self):
            raise TypeError('No __init__ is possible since no C++ constructors were declared')
    else:
        def __init__(self, *args, _new=init, return_type=None, signature=None):
            self.move_from(_new(*args, return_type=return_type, signature=signature))
    return __init__

################################################################################

def render_member(key, value, old):
    if old is None:
        def fget(self, _old=value):
            return _old(self)._set_ward(self)
    else:
        def fget(self, _old=value, _return=old):
            out = _old(self)._set_ward(self)
            return out.cast(_return)

    def fset(self, other, _old=value):
        _old(self).copy_from(other)

    name = getattr(old, '__name__', old)

    doc = 'C++ member variable' if name is None else 'Member of type `{}`'.format(name)
    return property(fget=fget, fset=fset, doc=doc)

################################################################################

def copy(self):
    '''Make a copy of self using the C++ copy constructor'''
    other = self.__new__(type(self))
    other.copy_from(self)
    return other

################################################################################

def find_class(mod, name):
    '''Find an already declared class and return its deduced properties'''
    old = common.unwrap(getattr(mod, name, None))
    annotations = {}
    props = {'__module__': mod.__name__}

    if old is not None:
        for c in reversed(old.mro()[:-1]): # skip base class 'object'
            props.update(c.__dict__)
            annotations.update(getattr(c, '__annotations__', {}))
    props['__annotations__'] = annotations
    return old, props

################################################################################

def render_type(translate, pkg: str, bases: tuple, name: str, methods):
    '''Define a new type in pkg'''
    mod, name = common.split_module(pkg, name)
    old_cls, props = find_class(mod, name)

    new = props.pop('__new__', None)
    if callable(new):
        log.warning('{}.{}.__new__ will not be rendered'.format(pkg, name))
    methods['__init__'] = render_init(methods.pop('new', None))

    for k, v in methods.items():
        k = common.translations.get(k, k)
        if k.startswith('.'):
            old = common.unwrap(props['__annotations__'].get(k[1:]))
            log.info("deriving member '%s.%s%s' from %s", mod.__name__, name, k, repr(old))
            translate['%s.%s' % (name, old)] = props[k[1:]] = render_member(k[1:], v, old)
        else:
            old = common.unwrap(props.get(k, common.default_methods.get(k)))
            log.info("deriving method '%s.%s.%s' from %s", mod.__name__, name, k, repr(old))
            translate[old] = props[k] = render_function(v, old)

    props.setdefault('copy', copy)

    translate[old_cls] = cls = type(name, bases, props)
    log.info("rendering class '%s.%s'", mod.__name__, name)
    setattr(mod, name, cls)
    return mod, cls

################################################################################

def render_object(pkg, key, value):
    '''put in the module level functions and objects'''
    mod, key = common.split_module(pkg, key)
    log.info("rendering object '%s.%s'", mod.__name__, key)
    localns = mod.__dict__

    old = common.unwrap(getattr(mod, key, None))
    if old is None:
        pass
    elif callable(value):
        log.info("deriving function '%s.%s' from %s", mod.__name__, key, repr(old))
        assert callable(old), 'expected annotation to be a function'
        value = render_function(value, old)
    else:
        log.info("deriving object '%s.%s' from %s", mod.__name__, key, repr(old))
        if not isinstance(old, type):
            raise TypeError('expected placeholder {} to be a type'.format(old))
        value = value.cast(old)

    setattr(mod, key, value)
    return old, value

################################################################################

def make_callback(origin, types):
    '''Make a callback that calls the wrapped one with arguments casted to the given types'''
    if not callable(origin):
        return origin
    def callback(*args):
        return origin(*(a.cast(t) for a, t in zip(args, types)))
    return callback

def is_callable_type(t):
    '''Detect whether a parameter to a C++ function is a callback'''
    t = getattr(t, '__origin__', None)
    return t is typing.Callable or t is collections.abc.Callable
# if not hasattr(typing, 'CallableMeta'):
#     raise ImportError('Python 3.7 has a bug where typing.CallableMeta is missing')

################################################################################

def render_function(fun, old):
    '''
    Replace a declared function's contents with the one from the document
    - Otherwise if the declared function takes 'function', call the declared function
    with the document function passed as 'function' and the other arguments filled appropiately
    - Otherwise, call the document function
    '''
    if old is None:
        def bound(*args, _orig=fun):
            return _orig(*args)
        return functools.update_wrapper(bound, common.opaque_signature)

    if isinstance(old, property):
        return property(render_function(fun, old.fget))

    sig = inspect.signature(old)

    has_fun = '_fun_' in sig.parameters
    if has_fun:
        sig = common.discard_parameter(sig, '_fun_')

    types = tuple(p.annotation.__args__ if is_callable_type(p.annotation) else None for p in sig.parameters.values())
    empty = inspect.Parameter.empty

    if '_old' in sig.parameters:
        raise ValueError('Function {} was already wrapped'.format(old))

    # Eventually all of the computation in wrap() could be moved into C++
    # (1) binding of args and kwargs with default arguments
    # (2) for each arg that is annotated with Callback, make a callback out of it
    # (3) either cast the return type or invoke the fun that is given
    # (4) ...
    if has_fun:
        def wrap(*args, _orig=fun, _bind=sig.bind, _old=old, gil=None, signature=None, **kwargs):
            bound = _bind(*args, **kwargs)
            bound.apply_defaults()
            # Convert args and kwargs separately
            args = (a if t is None or a is None else make_callback(a, t) for a, t in zip(bound.args, types))
            kwargs = {k: (v if t is None or v is None else make_callback(v, t)) for (k, v), t in zip(bound.kwargs.items(), types[len(bound.args):])}
            return _old(*args, _fun_=_orig, **kwargs)
    else:
        ret = sig.return_annotation

        for k, p in sig.parameters.items():
            if p.kind == p.VAR_KEYWORD or p.kind == p.VAR_POSITIONAL:
                raise TypeError('Parameter {} cannot be variadic (e.g. like *args or **kwargs)'.format(k))

        def wrap(*args, _orig=fun, _bind=sig.bind, _return=ret, gil=None, signature=None, **kwargs):
            bound = _bind(*args, **kwargs)
            bound.apply_defaults()
            # Convert any keyword arguments into positional arguments
            out = _orig(*(make_callback(a, t) if t is not None else a for a, t in zip(bound.arguments.values(), types)))
            if _return is empty:
                return out # no cast
            if _return is None or _return is type(None):
                return # return None regardless of output
            if out is None:
                raise TypeError('Expected {} but was returned object None'.format(_return))
            return out.cast(_return)

    return functools.update_wrapper(wrap, old)

