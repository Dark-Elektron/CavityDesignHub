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

Resonant cavities are particularly responsible for particle acceleration. Over the last
decade, there has been a great deal of advancement in the design and fabrication of
resonant cavities. Cavity shapes have evolved from being totally cylindrical in shape,
the so-called pillbox cavities, to the widely used (standard) elliptical cavities nowadays.
Cavity geometric parameters are selected based on certain desired properties. The so-
called cavity figures of merit are used to quantify these desired properties.

Figures of merit in accelerator cavity design
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are several figures of merits that quantify the efficiency and effectiveness of accelerating cavities.
These include the quality factor, :math:`Q`, the shunt impedance, :math:`R_\mathrm{sh}`, dissipated power,
:math:`P_\mathrm{c}`, stored energy, :math:`U`, etc.

The quality factor quantifies the ratio of the amount of energy stored in a cavity to the power dissipated at the
cavity walls. Mathematically,

.. math::
   Q = \omega \frac{U}{P_\mathrm{ds}}.

The shunt impedance is another important quantity used to characterise the losses in a cavity.
It is given mathematically as:

.. math::
   R_\mathrm{sh} = \frac{V_\mathrm{acc}^2}{P_\mathrm{ds}}

where :math:`U` is the stored energy given as:

.. math::
   U = \frac{1}{2} \epsilon_0 \int_V |E|^2 \mathrm{d}V = \frac{1}{2} \mu_0 \int_V |H|^2 \mathrm{d}V.

and :math:`P_\mathrm{c}` is the power loss given as:

.. math::
   P_\mathrm{ds} = \frac{1}{2} R_\mathrm{sh} \int_S |H|^2 \mathrm{d}s.

The accelerating voltage is given as

.. math::
   V_\mathrm{acc} = \int_0^L E_z (\rho, \phi) e^{j\omega t}\mathrm{d}z.

where L is the length of the cavity or beam pipe.


.. _QUICK:Electromagnetic field theory:

Electromagnetic field theory
*********************

Maxwell equations
Electromagnetic interactions are represented mathematically with the Maxwell equations.

.. math::

     \nabla D = \rho,

     \nabla B = 0,

     \nabla \times H = J + \frac{\partial D}{\partial t},

     \nabla \times E = - \frac{\partial B}{\partial t}.

where E is the electric field intensity, and H is the magnetic field intensity, D is the electric flux density,
B is the magnetic flux density. They are related by material properties given as follows:

.. math::

   D = \epsilon E,

   B = \mu H.

Combining the Maxwell equations and using appropriate vector identities result in the wave equations for electric
and magnetic fields as follows.

.. math::

    \left(\nabla^2 - \frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)\bigg\{\begin{matrix} E\\H\end{matrix}\bigg\} = 0.
magnetic permeability and :math:`\epsilon` electrical permittivity.


.. _QUICK:Maxwell eigenvalue problem:

Maxwell eigenvalue problem
^^^^^^^^^^^^^^^^^^^^

The Maxwell Eigenvalue Problem (MEVP) is solved using SLANS \cite{SLANS} to evaluate
:math:`e_\mathrm{pk}`, :math:`b_\mathrm{pk}`, and :math:`R/Q`. The MEVP is given as

.. math::
   \nabla \times \left({\mu}\, \nabla \times E(x)\right) - \lambda(x)\epsilon\,  E(x)= 0, & & \lambda = \frac{\omega^2}{c^2},~\mathbf{x} \in \mathbb{R}^7,

   \nabla \cdot E = 0 & & E \in \Omega,

   n \times E = 0 & & E \in \partial \Omega,


where:math:`\mathbf{x}` is a vector of the geometric variables describing the domain :math:`\Omega` with boundary
:math:`\partial \Omega`, :math:`\mathbf{E}` is the electric field, :math:`\lambda` a vector of eigenvalues, :math:`\mu`


.. _QUICK:Wakefield equations:

Wakefield analysis
^^^^^^^^^^^^^^^^

The longitudinal and transverse wake functions :math:`w_\parallel` and :math:`\mathbf{w}_\perp`, respectively,
are evaluated using ABCI electromagnetic code \cite{ABCI}. They are defined as

.. math::

   w_\parallel(\rho, s) &= -\frac{c}{q} \int E_z|_{z=ct-s} \mathrm{d}t,

   \mathbf{w}_\perp(\rho, s) &= \frac{c}{q} \int (\textbf{E}_\perp + c \hat z \times \textbf{B})|_{z=ct-s},

where :math:`s` is the distance between the leading and a trailing test particle with offset
:math:`\boldsymbol{\rho} = (x, y)` relative to the :math:`z`-axis, and :math:`z` is the direction of travel of the
particles. The longitudinal ($Z_\parallel$) and transverse impedances (:math:`Z_\perp`) are evaluated as the
Fourier transform of the wake functions thus:

.. math::
   Z_\parallel (\omega) &= \frac{1}{c} \int_0^\infty w_\parallel(s) \mathrm{e}^{(i\omega s/c)},

   \textbf{Z}_\perp (\omega) &= \frac{1}{c} \int_0^\infty \mathbf{w}_\perp(s) \mathrm{e}^{(i\omega s/c)}.



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
