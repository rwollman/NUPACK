

class ConversionError(TypeError):
    '''Default error type for all rebind type conversion errors'''


################################################################################

class Config:
    def __init__(self, methods):
        self._set_debug = methods['set_debug']
        self._get_debug = methods['debug']
        self.set_type_error = methods['set_type_error']
        self.set_type_names = methods['set_type_names']
        self.set_type = methods['set_type']
        self.set_output_conversion = methods['set_output_conversion']
        self.set_input_conversion = methods['set_input_conversion']
        self.set_translation = methods['set_translation']

    @property
    def debug(self):
        return self._get_debug().cast(bool)

    @debug.setter
    def debug(self, value):
        self._set_debug(bool(value))

################################################################################

from .render import render_module, render_init, render_member, \
                    render_type, render_object, render_function

from .common import update_module

from .types import Tuple, List, Dict, Callable, Array, Union

################################################################################

_declared_objects = []

def forward(x):
    _declared_objects.append(x)
    return x