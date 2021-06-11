# Flambaum-Ginges radiative potential

A c++ [c++-17] implementation of the Flambaum-Ginges radiative potential method, used for including quantum electrodynamics corrections into calculations for many-electron atoms.

First introduced in:
 * V. V. Flambaum and J. S. M. Ginges, Phys. Rev. A 72, 052115 (2005).

With finite-nuclear-size corrections, and updated fitting factors, from;
 * J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
 * J. S. M. Ginges and J. C. Berengut, J. Phys. B 49, 095001 (2016).

Note:
 * All radii (r and rN) are given in atomic units
 * Sign convention: H -> H - V_rad
 * High- and low- frequency V(r) do not include A or B fitting coefficients.
 * V_rad = V_Uehling + A(Z,l)*V_SEh + B(Z,l)*V_SEl + V_WK + i(g.n)H_Magnetic

See provided FGRP.pdf for brief definitions; see above papers for full explanation.

## Usage:

 * Just add the FGRadPot.cpp and FGRadPot.hpp files to your code source directly
 * See the example program for an example
 * Alternatively, use the example program to generate and output radiative potential
 * Requires GSL (GNU scientific libraries) https://www.gnu.org/software/gsl/
   * Install GSL libraries on ubuntu: _$sudo apt install libgsl-dev_
   * Install GSL libraries on mac (homebrew): _$brew install gsl_
 * On mac, may have to explicitely link to GSL; typically installed here:
   * _-I/usr/local/opt/gnu-scientific-library/include/_
   * _-L/usr/local/opt/gnu-scientific-library/lib/_
 * Uses c++-17 (g++ version 7 or higher; clang++ version 6 or higher)


## Functions:

```cpp
//! Uehling potential (r, rN in atomic units)
double V_Uehling(double Z, double r, double rN = 0.0);

//! Magnetic form-factor (r, rN in atomic units)
double H_Magnetic(double Z, double r, double rN = 0.0);

//! High-freq electric SE (NOT including Al) (r, rN in atomic units)
double V_SEh(double Z, double r, double rN = 0.0);

//! Low-freq electric SE (NOT including Bl) (r, rN in atomic units)
double V_SEl(double Z, double r, double rN = 0.0);

//! Effective Wickman-Kroll; not including FNS corrections
double V_WK(double Z, double r, double);

//! Al(Z) fitting function [PRA 93, 052509 (2016)]
double Fit::A(double Z, int l = 0);

//! Bl(Z) fitting function [PRA 93, 052509 (2016)]
double Fit::B(double Z, int l = 0);

```

 * Z is nuclear charge number
 * r is radial coordinate, in atomic units (aB=1)
 * rN is nuclear radius, in atomic units
   * nb: rN = Sqrt(5/3)*r_rms  [r_rms is root-mean-square radius]
