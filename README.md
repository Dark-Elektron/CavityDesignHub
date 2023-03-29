![GitHub all releases](https://img.shields.io/github/downloads/Dark-Elektron/CavityDesignHub/total?logo=Github) 
![GitHub issues](https://img.shields.io/github/issues-raw/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub pull requests](https://img.shields.io/github/issues-pr/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub closed pull requests](https://img.shields.io/github/issues-pr-closed-raw/Dark-Elektron/CavityDesignHub?logo=Github)


Overview
=======

This repository contains Python codes for conducting analysis on accelerating
cavities. Eigenmode analysis, wakefield analysis, and general post-processing.

Each module performs a different operation. The analysis that are currently
supported in this module are eigenmode analysis, wakefield analysis,
and multipacting analysis.

* Eigenmode analysis - SLANS [[1]](#1)
* Wakefield analysis - ABCI [[2]](2#)
* Optimisation - Python
* Uncertainty quantification - Python
* Postprocessing - Python

## Eigenmode Analysis

Eigenmode analysis is performed using the SLANS electromagnetic code. The code
also calculates most of the figures of merit. Some postprocessing is, however,
required to transform them to the form that is used in most papers related
to accelerating cavities design.

The SUPERLANS code is intended to calculate azimuthal-homogenous modes in
axissymmetric cavities, periodical structure, and cut-off frequencies in
long homogenous waveguides. SLANS is written by Sergey
Belomestnykh and it consists of a set of executable files for different
purposes.


## Wakefield Analysis

Wakefield analysis is performed using the ABCI electromagnetic code which solves the Maxwell
equations directly in the time domain when a bunched beam goes through an axi-symmetric
structure on or off axis. An arbitrary charge distribution can be defined 
by the user (default=Gaussian)



## References
   <a id="1">[1]<a/>
   D. Myakishev and V. Yakovlev, "The new possibilities of SUPERLANS code for evaluation of 
   axisymmetric cavities", in Proc. of the 1995 Particle Accelerator Conf. 
   Dallas, TX, 1995, pp. 2348-50, Available: 
   http://epaper.kek.jp/p95/ARTICLES/MPC/MPC17.PDF.

   <a id="2">[2]<a/>
   Y. H. Chin, Azimuthal Beam Cavity Interaction (ABCI), https://abci.kek.jp/
   