from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import pybind11


class get_pybind_include(object):
    """Helper to get pybind11 include path"""
    def __str__(self):
        return pybind11.get_include()

ext_modules = [
    Extension(
        name="sph_cpp",  # The name of the Python module
        sources=["sph_bindings.cpp"],  # Add more .cpp files if needed
        include_dirs=[
            get_pybind_include(),
            ".",  # current directory if sph.cpp is here
        ],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17"],  # Or c++14 if needed
    )
]

setup(
    name="sph_cpp",
    version="0.1.0",
    author="Your Name",
    description="SPH Simulation Python Binding",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
