![GitHub all releases](https://img.shields.io/github/downloads/Dark-Elektron/CavityDesignHub/total?logo=Github) 
![GitHub issues](https://img.shields.io/github/issues-raw/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub pull requests](https://img.shields.io/github/issues-pr/Dark-Elektron/CavityDesignHub?logo=Github) 
![GitHub closed pull requests](https://img.shields.io/github/issues-pr-closed-raw/Dark-Elektron/CavityDesignHub?logo=Github)


Overview
=======

This is the introduction file and I have to write somethign at some point in time here
Now, this software is used to for conducting analysis on accelerating
cavities. Eigenmode analysis, wakefield analysis, multipacting analysis
and general post-processing.

Each module performs a different operation. The analysis that are currently
supported in this module are eigenmode analysis, wakefield analysis,
and multipacting analysis.

* Eigenmode analysis - SLANS :cite:p:`SLANS`
* Wakefield analysis - ABCI :cite:p:`ABCI`
* Multipacting analysis - Multipac
* Optimisation - Python
* Uncertainty quantification - Python
* Postprocessing - Python

Eigenmode Analysis
***************

Eigenmode analysis is performed using the SLANS electromagnetic code. The code
also calculates most of the figures of merit. Some postprocessing is, however,
required to transform them to the form that is used in most papers related
to accelerating cavities design.

The SUPERLANS code is intended to calculate azimuthal-homogenous modes in
axissymmetric cavities, periodical structure, and cut-off frequencies in
long homogenous waveguides :cite:p:`SLANS`. SLANS is written by Sergey
Belomestnykh and it consists of a set of executable files for differnt
purposes. The first of these is the ``genmesh.exe`` which reads a geometry
file ``<filename>.geo`` written using Python and generates the mesh file
which is a .gem file and some other related files. ``slansre.exe`` is then
called to run the eigenmode simulation and the results are output to specified folder.

The files output by the SLANS codes are basically three types

* binary files:
* text files: which can be read by regular text editors
* meta files

The inputs and output files of the various SLANS codes are given below:



Wakefield Analysis
******************



.. bibliography::
   :all: