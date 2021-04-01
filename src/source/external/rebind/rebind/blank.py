

def _blank():
    pass

class _Blank:
    def _blank(self):
        pass

_blank_function = _blank.__code__.co_code
_blank_method = _Blank._blank.__code__.co_code

def check_rendering(name, obj):
    # types.FunctionType
    if isinstance(obj, types.FunctionType) and hasattr(obj, '__code__'):
        if obj.__code__.co_code == _blank_function:
            raise ImportError(name + ' is an empty function')
        if obj.__code__.co_code == _blank_method:
            raise ImportError(name + ' is an empty method')

    if isinstance(obj, types.FunctionType) and hasattr(obj, '__dict__'):
        for k, v in obj.__dict__.items():
            check_rendering(name + '.' + k, v)

