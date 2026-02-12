# hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System

<!-- ![hi_class logo](docs/hi_class_logo.gif) -->
<!-- <img src="docs/hi_class_logo.gif" alt="hi_class logo" width="140" align="right"> -->
<!-- <img src="docs/hi_class_logo.gif" alt="hi_class logo" width="140"> -->

hi_class extends the CLASS Boltzmann code to cover Horndeski and related scalar-tensor models of dark energy and modified gravity. It is based on CLASS by Julien Lesgourgues, with major inputs from Thomas Tram and others.

- Website: http://hiclass-code.net
- CLASS website: http://class-code.net

## Authors

- Miguel Zumalacarregui
- Emilio Bellini
- Ignacy Sawicki

## Installation

### From PyPI/TestPyPI

```bash
pip install hiclassy
```

This installs the Python wrapper and builds the C core. You need a working C compiler (e.g. gcc) and, optionally, OpenMP support for parallel execution.

### From source

```bash
make clean
make class
```

To build the Python wrapper, run:

```bash
make
```

If compilation fails, check the Makefile for compiler, optimization flags, and OpenMP settings.

## Quick start

Use the Python interface:

```python
from hiclassy import HiClass
```

The HiClass object is the equivalent of the Class object for standard Class. All methods and attributes, plus additional HiClass specific, are shared between the two. Then, the usage of the two should be equivalent. 

If you want to use the C executable instead, you can run:

```bash
./class explanatory.ini
```

Parameter documentation and examples are available in `hi_class.ini` and `explanatory.ini`, plus the example files in `gravity_models/`.

## Citing hi_class

If you use hi_class, please cite:

- M. Zumalacarregui, E. Bellini, I. Sawicki, J. Lesgourgues, P. Ferreira, "hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System", JCAP 1708 (2017) no.08, 019, http://arxiv.org/abs/arXiv:1605.06102
- E. Bellini, I. Sawicki, M. Zumalacarregui, "hi_class: Background Evolution, Initial Conditions and Approximation Schemes", http://arxiv.org/abs/arXiv:1909.01828

Please also cite the relevant CLASS papers, including:

- CLASS I: Overview, http://arxiv.org/abs/1104.2932
- CLASS II: Approximation schemes, http://arxiv.org/abs/1104.2933

## Plotting utilities

The package includes the Class Plotting Utility `CPU.py` for plotting Cl's, P(k), and related outputs, including model comparisons. Run:

```bash
python CPU.py --help
```

A MATLAB helper is available in `plot_CLASS_output.m`.

## Development

We recommend developing from the GitHub repository:

https://github.com/emiliobellini/hi_class_public

For hi_class-specific updates, see this repository and the `gravity_models/` examples.

## Support

For support, please open an issue in the repository or refer to the documentation at http://hiclass-code.net.
