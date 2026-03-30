# myPackage/__init__.py
import pyomo.environ as pe

# Importing from the main package modules
from .file_management import *

# from .special_add_loopless_implementation import *
# from .special_fastcc_implementation import *

# Importing the Legacy sub-package
from .Legacy import *

# Optional: Define an __all__ variable to specify what is available for import when using 'from myPackage import *'
__all__ = file_management.__all__ + Legacy.__all__
