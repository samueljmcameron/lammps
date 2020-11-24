.. index:: compute active pressure

compute active pressure command
========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID activepressure keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* activepressure = style name of this compute command
* zero or more keywords may be appended
* keyword = *wrapcoords*

Examples
""""""""


.. parsed-literal::

   compute 1 all activepressure
   compute 1 all pressure NULL pair bond
   compute 1 all pressure NULL pair/hybrid lj/cut

Description
"""""""""""

Define a computation that calculates the active pressure of the entire system
of atoms.  The specified group must be "all". For now, this compute is only
valid for Active Brownian Particle systems in periodic boundary conditions,
i.e. systems which have Langevin equations of motion

.. math::

   \frac{d r_i}{dt} = v_0 e_i + \frac{F_i}{\gamma} + \sqrt{2D^{(t)}}\xi_i^{(t)}

where v_0 is the self propulsion velocity, F_i is some conservative force,
gamma is a drag coefficient (real units of kg/s), D^{(t)} is the bare (i.e.
with v\_0 = 0) translational diffusion coefficient, and xi is a noise with
unit variance and zero mean (real units of s^{-1/2}). The orientation of
the spheres, denoted e\_i, obey a rotational diffusion equation with diffusion
coefficient D^{(r)} (units of 1/s). 

The pressure is computed by the formula

.. math::

   P = \frac{\gamma}{d V}\sum_i^Nr_i \bullet (\sqrt{2D^{(t)}}\xi_i^{(t)} + v_0 e_i )

and so only computes the active swim pressure (to get full pressure, you must
additionally compute pressure as usual but neglect kinetic energy). N is the
number of atoms in the system,  d is the
dimensionality of the system (2 or 3 for 2d/3d), and V is the system volume
(or area in 2d).  When periodic boundary conditions are used, it is not clear
whether the unwrapped coordinates vs image coordinates for r_i should be used
in the equation above, so both will be tested. The above equation is based off
of work from ref. :ref:`(Winkler) <Winkler1>`, but does not assume that the
system is in steady state (compare to e.g. eqn (44) of ref.
:ref:`(Winkler) <Winkler1>`, which assumes steady state via eqn (31) of the
same paper).

If no extra keyword is listed, the entire equation above is
calculated using unwrapped coordinates for the position vectors r_i.
If *wrapcoords* is listed, then the position vectors r_i will be
in their wrapped state (i.e. r_i will be within the periodic box
boundaries).

No temperature is included in this calculation, as temperature is not a
well defined quantity for active systems. Instead, the pressure due to
translational diffusion is included in the first term of the pressure
formula above.


----------


**Output info:**

This compute calculates a global scalar (the active pressure).
These values can be used by any command that uses global scalar
or vector values from a compute as input.  See the
:doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS
output options.

The scalar value calculated by this compute is
"intensive".  The scalar value will be in pressure
:doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute pressure <compute_pressure>`,

**Default:** none


----------

.. _Winkler1:



**(Winkler)** Winkler, Wysocki, Gompper, Soft Matter, 11, 6680 (2015).
