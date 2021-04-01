from .cpp import document
from rebind.build import render_module
from . import submodule

rendered_document = render_module(__name__, document)
