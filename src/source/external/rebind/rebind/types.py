from typing import Tuple, List, Dict, Callable, Union, Optional

class ArrayType:
    def __getitem__(self, type):
        assert not isinstance(type, tuple), 'Only 1 argument may be given to Array[]'
        return Tuple[type, ...]

Array = ArrayType()

class Forward:
    def __init__(self, placeholder):
        self._origin_ = placeholder

    def __repr__(self):
        return 'Forward({!r})'.format(self._origin_)

    def __getattr__(self, key):
        raise TypeError('Undefined attribute {!r} of forward object {!r}'.format(key, self._origin_))

    def __call__(self, *args, **kws):
        raise TypeError('Forward object {!r} cannot be called'.format(self._origin_))

# also: None, bool, int, float, str, bytes, object

# alternative syntax?: e.g. (T, U) instead of Tuple[T, U]
# Any? Iterator? NamedTuple?
