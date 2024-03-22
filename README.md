
# SCARCCpy - Simulation of Combined Antibiotics in Cross-feeding Communities

SCARCC is a python package developed by the Harcombe Lab designed to identify synergistic interactions between antibiotics within our two-species cross-feeding microbial community. 

This package utilizes [COBRA](https://github.com/opencobra/cobrapy) to analyze the perturbation effects caused by antibiotics in the genome-scale metabolic network using Flux Balance Analysis (FBA), and [COMETS](https://github.com/segrelab/cometspy) simulations to incorporate multispecies growth simulation using Dynamic Flux Balance Analysis (dFBA). 


# Documentation
The documentation is at [readthedocs](https://scarccpy.readthedocs.io/en/latest/index.html)


# Installation

[Python version](https://www.python.org/downloads/release/python-3100/) greater than Python 3.10 is required.

Use pip to install, use of [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/) is encouraged. 

For more detail setting up in MSI, see [INSTALL.rst](https://github.com/TFwongw/scarccpy/blob/main/INSTALL.rst)

```py
pip install scarcc
```
# License
The SCARCC source code is available under the MIT License.