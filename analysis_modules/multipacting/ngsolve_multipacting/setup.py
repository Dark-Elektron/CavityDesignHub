from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        "hit_bound_cy.pyx", compiler_directives={"language_level": "3"}, annotate=True
    )
)