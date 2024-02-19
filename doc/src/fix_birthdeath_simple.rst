.. index:: fix birthdeath/simple
.. index:: fix birthdeath/simple/1D/ratchet
.. index:: fix birthdeath/simple/1D/smoothratchet

fix birthdeath/simple command
=============================

fix birthdeath/simple/1D/ratchet command
========================================

fix birthdeath/simple/1D/smoothratchet command
==============================================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style_name alive_type dead_type seed birthrate deathrate shift cleanevery keyword args

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *birthdeath/simple* or *birthdeath/simple/1D/ratchet* or *birthdeath/simple/1D/smoothratchet*
* alive_type = type of atom which is considered alive
* dead_type = type of atom which is considered dead
* birthrate = rate of atom creation 
* deathrate = rate of atom annihilation
* shift = the total distance between the two daughter atoms
* cleanevery = how often to delete dead_type atoms (does not affect results, just saves space)
* one or more keyword/value pairs may be appended
* keyword (must be in order if included) = *exactpoisson* or *omega* or *fraction* or *nterms* or *height*

  .. parsed-literal::

     *exactpoisson* value = *yes* or *no*
       *yes* = a maximum of one event (birth or death) can happen per time step
       *no* = each particle could give birth or die at each time step (default)
     *omega* value = *omega*
       *omega* = geometrical projection factor to scale the birth projection
     *fraction* value = *fraction*
       *fraction* = fraction of the 1D domain where the ratchet potential is decreasing
     *nterms* value = *nterms*
       *nterms* = number of terms to include in the smooth ratchet potential series (more terms means more asymmetry)
     *height* value = *height*
       *height* = the height of the ratchet potential (in units of energy)

       
Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all birthdeath/simple 1 2 12908410 0.04 0.0001 0.3 500
   fix 1 all birthdeath/simple 1 2 12908410 0.04 0.0001 0.3 500 exactpoisson yes
   fix 1 all birthdeath/simple/1D/ratchet 1 2 12908410 0.04 0.0001 0.3 500 exactpoisson yes omega 1.0 fraction 0.25 height 1.0
   fix 1 all birthdeath/simple/1D/ratchet 1 2 12908410 0.04 0.0001 0.3 500 omega 1.0 fraction 0.25 height 1.0
   fix 1 all birthdeath/simple/1D/smoothratchet 1 2 12908410 0.04 0.0001 0.3 500 omega 1.0 nterms 10 height 1.0
   fix 1 all birthdeath/simple/1D/smoothratchet 1 2 12908410 0.04 0.0001 0.3 500 exactpoisson yes omega 1.0 nterms 10 height 1.0

   
Description
"""""""""""


Allow the atoms to give birth or die at each time step. If *exactpoisson* = *yes*, then this
means that a single atom (randomly selected at each timestep) can either make a copy of
itself with probability *birthrate x number of particles x dt*, die with probability
*deathrate x number of particles x (number of particles -1) x dt*, or do neither).
If *exactpoisson* = *no* (default), then each atom can make a copy of itself with probability
*birthrate x dt*, die with probability *deathrate x (number of particles -1) x dt*, or
do neither. Either process ensures that the total number of particles in steady-state
will be *birthrate / deathrate*. This birth-death process is implemented in the
``post_integrate()`` method (see how a timestep works for details of when this method is called
in a timestep).

If *exactpoisson* = *yes*, then the population dynamics of the system (once integrating out
spatial degrees of freedom) is described by the master equation

.. math::

   \frac{ d p_n(t)}{ dt} = k_0 n(n+1)p_{n+1}(t) - k_0n(n-1)p_n(t) + b(n-1)p_{n-1}(t) - bnp_{n}(t)

where :math:`b` = *birthrate*, :math:`k_0` = *deathrate*, :math:`n` = *number of particles*
and :math:`p_n(t)` is the probability of :math:`n` particles existing at time :math:`t`. If
*exactpoisson* = *no*, the process is more difficult since multiple events may occur in
a given timestep.

When an atom makes a copy of itself, the positions of the two daughter atoms (the original and
the copy) depend on which fix style is being used. For *birthdeath/simple*, the daughter
atoms maintain the same centroid as prior to their split, but are separated by a distance
*shift* along a randomly oriented vector in space. For *birthdeath/1D* styles, the
daughter atoms are shifted via the equation

.. math::

   x_{new} = x_{old} \pm \frac{shift}{2}(1+(\partial_x U)^2/\omega^2)^{-1/2}

where :math:`-\partial_x U` is the external force (a piecewise ratchet potential
for *birthdeath/1D/ratchet* and a smooth ratchet potential for *birthdeath/1D/smoothratchet*)
evaluated at :math:`x_{old}` (which is the position of the parent atom before it splits).

In all cases, the velocities of the daughter atoms are half that of the parent atom.


For the *birthdeath/1D* styles, an external force :math:`-\partial_x U` is also applied to all
the particles (a piecewise ratchet potential for *birthdeath/1D/ratchet* and a
smooth ratchet potential for *birthdeath/1D/smoothratchet*). The force is applied using the
``post_force()`` method (see how a timestep works for details of when this method is called
in a timestep).





---------

.. note::

   These fixes are likely best used with the overdamped :doc:`fix brownian <fix_brownian>`
   dynamics using the ``final_integrate`` option (so that births happen prior to
   time integration).
   

----------


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BIRTHDEATH package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix brownian <fix_brownian>`

Default
"""""""

The default for *exactpoisson* is *no*.
