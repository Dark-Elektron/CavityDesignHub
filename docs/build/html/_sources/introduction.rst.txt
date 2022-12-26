########
Overview
########

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


.. Note::

    It is planned that in future releases, all analysis codes will be custom
    codes written in Python.

Wakefield Analysis
******************

Wakefield analysis is performed using the ABCI
(Azimuthal Beam Cavity Interaction) code written by Yang Ho Chin
:cite:p:`ABCI`. It solves the Maxwell equations directly in the tiem domain
when a bunched beam goes through an axi-symmetric structure on or off axis.
It can be found `here <https://abci.kek.jp/abci.htm>`_.

.. Note::

    It is planned that in future releases, all analysis codes will be custom
    codes written in Python.

Multipacting Analysis
********************

Multipacting simulations are performed using Multipac code by . The software
can be obtained as described
`here <https://accelconf.web.cern.ch/e08/papers/mopp137.pdf>`_.
Multipac is a set of executable files and MATLAB codes.

.. Note::

    It is planned that in future releases, all analysis codes will be custom
    codes written in Python.

Optimisation
***********

Optimisation is done using self-written Python codes which call to the previously
mentioned analysis codes. Currently, brute force and Genetic Algorithm (GA)
optimisations are supported. More should be included in future releases.

Uncertainty Quantification
***********

Uncertainty quantificaion (UQ) is done using self-written Python codes
which call to the previously mentioned analysis codes. The procedure followed
is described in here. The mathematics can be found in the theory section.

.. bibliography::
   :all: