from setuptools import find_packages, setup
from pathlib import Path

setup(
    name='pde-rk',
    version='0.1.0',
    license="CC BY 4.0",
    author='Tom Bland',
    author_email='tom_bland@hotmail.co.uk',
    packages=find_packages(),
    install_requires=['numpy',
                      'matplotlib',
                      'scipy',
                      'pandas'],
    description='Package for simulating one-dimensional PDE models with and adaptive Runge-Kutta scheme'
)
