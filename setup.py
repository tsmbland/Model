from setuptools import find_packages, setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pde-rk",
    version="0.1.4",
    license="CC BY 4.0",
    author="Tom Bland",
    author_email="tom_bland@hotmail.co.uk",
    packages=find_packages(),
    install_requires=["numpy"],
    description="Package for simulating one-dimensional PDE models with and adaptive Runge-Kutta scheme",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Source Code": "https://github.com/tsmbland/pde-rk",
    },
)
