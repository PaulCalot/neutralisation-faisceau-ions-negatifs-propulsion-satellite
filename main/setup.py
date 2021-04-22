from distutils.core import setup
from Cython.Build import cythonize
import numpy
setup(
    ext_modules=cythonize("handler_collision_with_wall.pyx"), # "DSMC.py"), # "DSMC.py"
    include_dirs=[numpy.get_include()]
)
