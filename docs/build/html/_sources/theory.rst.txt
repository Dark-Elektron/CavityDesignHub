########
Theory
########

In this section the theory governing the various analysis in this program are
introduces. Note that this would not be exhaustive as that would require several
semester's worth of instruction. What is attempted here is to give the fundamental
equations, some basic descriptions and method in which the code was written so that
one can understand the different line of codes. This is done so that the user can
better understand and interpret exactly what results the program outputs and perhaps
point out inconsistencies in the theory and code or even ways in which a routine can
be optimised for speed and accuracy.

.. contents:: Contents of this Page
   :local:


.. _QUICK:Accelerator cavity basics:

Aceclerator cavity basics
**************************
Particle accelerators are devices used for the acceleration of charged particles to high
velocities to serve for various purposes. This is achieved using electromagnetic fields.
They are specially used for fundamental science research in areas cutting across physics,
medicine, material science, and so on. Typical applications today include high energy
physics research, biomedical research, radiotherapy, and ion implantation for the manu-
facturing of semiconductors. Particle accelerators are comprised of several components
which include the particle or beam source, a radiofrequency (RF) source, electromagnets,
waveguides, resonant cavities, etc.

.. _ellipse tangent:

.. figure:: ../images/ellipse_tangent.png
   :alt: accelerator cavity
   :align: center
   :width: 200px

Resonant cavities are particularly responsible for particle acceleration. Over the last
decade, there has been a great deal of advancement in the design and fabrication of
resonant cavities. Cavity shapes have evolved from being totally cylindrical in shape,
the so-called pillbox cavities, to the widely used (standard) elliptical cavities nowadays.
Cavity geometric parameters are selected based on certain desired properties. The so-
called cavity figures of merit are used to quantify these desired properties.

.. _QUICK:Electromagnetic field theory:

Electromagnetic field theory
*********************


.. _QUICK:Maxwell eigenvalue problem:

Maxwell eigenvalue problem
*********************


.. _QUICK:Wakefield equations:

Wakefield analysis
*********************



.. _QUICK:Theory of multipacting:

Theory of multipacting
**********************

Multipacting is a phenomenon in resonant frequency structures whereby charged
particles are continuously discharged at an exponential rate from the
conductor walls of the device. Multipacting occurs only at specific
conditions dependent on the alternating field and the wall's surface
properties. More specifically, the emission properties of the wall coupled
with the electromagnetic field strength and profile are determinants of
multipacting. In accelerating cavities, multipacting occurs mostly at the
equator.

For multipacting to occur, two conditions must be satisfied. The first is
that the an electron emitted at the cavity wall and driven by the cavity
electromagnetic field returns to the same point after an integer number of
cycles. Secondly, the impacting electron produces more than one secondary
electron \cite{Oija}. Different measures could predict Multipacting. The one
adopted here is implemented in the 2D code described in \cite{Oija, Padamse}.
The tools for multipacting analysis are the counter and distance functions.
The counter and distance functions are used to predict multipacting.

The counter functions are the electron counter function :math:`c_N`, enhanced
counter function :math:`e_N` and total electron counter function :math:`t_N`. These
quantify the number of free electrons remaining after a given number of
impacts (20 in this case), the number of secondary electrons and the number
of all electrons, respectively, where :math:`N` is the maximum number of impacts.
The final impact energy :math:`Ef_N` of the free electrons is also calculated.

The distance function :math:`d_N` defined mathematically as

.. math::

    d_N = \sqrt{|x_N-x_0|^2 + \gamma^2|\mathrm{e}^{i\psi_N} - \mathrm{e}^{i\psi_0}|^2}

quantifies the distance between the initial point and the last impact point where
:math:`(x_0, \psi_0)`, :math:`(x_N, \psi_N)` are the coordinates of the initial and
final position after :math:`N` impacts.

.. figure:: ../images/sey_Nb.png
   :alt: accelerator cavity
   :align: center
   :width: 200px

Multipacting strongly depends on the surface finish of the cavity wall
materials and, therefore, strongly depends on the material's secondary
emission yield curve (SEY). It is quite peculiar because if the energy of
the impacting particle is too low or too high, no electron is released.
This sets an experimentally determined bound of impact energy for which
electrons are released in different materials. Figure \ref{fig: sey niobium}
shows a typical SEY curve for Niobium.


.. _QUICK:Sensitivity analysis and UQ:

Sensitivity analysis and uncertainty quantification
*********************


.. _QUICK:Cavity optimization:

Elliptical Cavity optimization
*********************

Accelerator cavities consists of cells joined together at either the equator
or the iris. Figure {} shows a typical elliptical accelerator cavity.
The geometry could be divided into two groups: the mid-cells and the end-cells
group. For optimisation purposes, it is computationally cheaper to optimise the
groups independently than optimising the entire cavity geometry. Figure {} shows
a typical parametrisation of an accelerator cavity. It is important to note that
several designers may parametrise the cavity in different ways ref{} and also
use different notations.
