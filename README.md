# hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System

<!-- ![hi_class logo](docs/hi_class_logo.gif) -->
<!-- <img src="docs/hi_class_logo.gif" alt="hi_class logo" width="140" align="right"> -->
<!-- <img src="docs/hi_class_logo.gif" alt="hi_class logo" width="140"> -->

hi_class extends the CLASS Boltzmann code to cover Horndeski and related scalar-tensor models of dark energy and modified gravity. It is based on CLASS by Julien Lesgourgues, with major inputs from Thomas Tram and others.

hi_class retains the functionality of the corresponding CLASS version and adds support for modified gravity. hi_class v3.4.0.0 is based on CLASS v3.4.0, so the standard CLASS v3.4.0 features, input parameters and outputs are available.

- Website: https://hiclass-code.net
- Documentation: https://github.com/hiclass-code/hi_class_public/wiki
- CLASS website: http://class-code.net

## Authors

- Emilio Bellini
- Ignacy Sawicki
- Miguel Zumalacarregui

## Installation

### From PyPI

```bash
pip install hiclassy
```

This installs the Python wrapper and builds the C core. You need Python 3.9 or newer, `make`, and C and C++ compilers (e.g. `gcc` and `g++`). Pip installs the Python dependencies and compiles the source package locally.

### From source

```bash
git clone https://github.com/hiclass-code/hi_class_public.git
cd hi_class_public
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

Since 9 February 2026, the wrapper uses `hiclassy` and `HiClass`. To adapt a standard CLASS Python script, replace `from classy import Class` with `from hiclassy import HiClass` and `Class()` with `HiClass()`. The standard methods of the corresponding CLASS version keep the same names and calling conventions; `HiClass` adds methods for modified gravity. See the [Python wrapper guide](https://github.com/hiclass-code/hi_class_public/wiki/Python-wrapper) for examples and conventions.

If you want to use the C executable instead, you can run:

```bash
./class explanatory.ini
```

Parameter documentation and examples are available in `hi_class.ini` and `explanatory.ini`, plus the example files in `gravity_models/`.

## Citing hi_class

If you use hi_class, please cite:

- M. Zumalacarregui, E. Bellini, I. Sawicki, J. Lesgourgues, P. Ferreira, "hi_class: Horndeski in the Cosmic Linear Anisotropy Solving System", JCAP 1708 (2017) no.08, 019, https://arxiv.org/abs/1605.06102
- E. Bellini, I. Sawicki, M. Zumalacarregui, "hi_class: Background Evolution, Initial Conditions and Approximation Schemes", https://arxiv.org/abs/1909.01828

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

https://github.com/hiclass-code/hi_class_public

For hi_class-specific updates, see this repository and the `gravity_models/` examples.

## Support

For support, please open an issue in the repository or refer to the [wiki](https://github.com/hiclass-code/hi_class_public/wiki).

## Weyl power and growth conventions

Write `W = (phi + psi)/2` and `Q_W(k,z) = k**4 P_W(k,z)`.
HiClass returns `Q_W` in `1/Mpc`, with k in `1/Mpc`. It is the power
spectrum of `k**2 W`, not the dimensionless power per logarithmic interval.

The clients `emu_like` and `hi_fast` accept q in `h/Mpc` and retain the
historical numerical normalization `S_W(q,z) = h**3 Q_W(h*q,z)` for
compatibility with existing datasets. This is not matter power in
`(Mpc/h)**3`; new `emu_like` Weyl headers label it `h^3/Mpc`.
Ratio targets additionally divide by their stored reference spectrum.

For every species, `f = (1/2) d ln P / d ln a = -(1+z)/(2P) dP/dz`.
For Weyl, P denotes the rescaled Weyl power; its growth can be negative.
Both clients differentiate HiClass power with second-order differences and
`dz = 1e-3`. They use a forward stencil near z=0. `emu_like` also uses a
backward stencil at the native upper time boundary; `hi_fast` reserves
coverage for its upper stencil. The fixed h normalization cancels in f.

Weyl requires a HiClass build exposing `pk_weyl`, `pk_weyl_lin`,
`get_pk_weyl`, and `get_pk_weyl_lin`. Both clients request `wPk`
automatically. `emu_like` follows the configured nonlinear setting and
rejects nonlinear Weyl requests; `hi_fast` explicitly uses linear power.
The unchanged `get_Weyl_pk_and_k_and_z` remains a comparison accessor, not
the source of client-side extrapolation.

Below the native k_min, the leading-order prescription is
`Q_W(k,z) = (k/k_min)**n_s Q_W(k_min,z)`. It assumes a k-independent Weyl
source at leading order and a single adiabatic analytic primordial power
law with `alpha_s = beta_s = 0`. Other primordial setups remain usable
within the native grid but are rejected below k_min. There is no high-k
or redshift extrapolation, and k must be strictly positive. Agreement
between implementations does not by itself establish this asymptotic
approximation's accuracy for every modified-gravity model.

Regenerating datasets updates targets and reference tables; existing
trained emulator weights are not changed by these code updates.
