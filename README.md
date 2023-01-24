# pde-rk

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsmbland/pde-rk/HEAD?filepath=%2Fscripts/simulate_par.ipynb)
[![CC BY 4.0][cc-by-shield]][cc-by]
[![PyPi version](https://badgen.net/pypi/v/pde-rk/)](https://pypi.org/project/pde-rk)

Functions for simulating PDE models with an adaptive Runge-Kutta scheme

## Installation

    pip install pde-rk

## Instructions

The repository contains an example notebook that demonstrates the method using the PAR polarity model (Goehring et al., 2011)

To run in the cloud, click 'launch binder' above.

To run on your local machine, follow these steps:

&#8291;1. Clone the repository:

    git clone https://github.com/tsmbland/pde-rk.git
    cd pde-rk

&#8291;2. Create conda environment:

    conda env create -f environment.yml

&#8291;3. Activate conda environment:

    conda activate pde-rk

&#8291;4. Open jupyter notebook:

    jupyter notebook scripts/simulate_par.ipynb

## License

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/

[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png

[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

