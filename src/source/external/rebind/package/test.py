def g(x,y,z): pass


def f(*args):
   print(*args)


import functools

sig = inspect.signature(g)

def run(*args, _orig=f, _bind=sig.bind, **kwargs):
    return _orig(*_bind(*args, **kwargs).args)

functools.update_wrapper(run, g)