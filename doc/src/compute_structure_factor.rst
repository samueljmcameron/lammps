.. index:: compute structurefactor

compute structure/sincos command
===================

compute structure/sincos/single command
===================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID style bin

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *structure/sincos* or *structure/sincos/single* = style name of this compute command
* bin = number of Fourier space (:math:`q`) bins (for *structure/sincos*) or which bin to calculate at (for *structure/sincos/single*)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all structure/sincos 100
   compute 1 all structure/sincos/single 3
   

Description
"""""""""""

Define a computation that calculates four quantities which are all
related to the spherically symmetric intermediate scattering function

.. math::

   F(q,t) = \left\langle\frac{1}{N}\sum_{j=1}^N\sum_{k=1}^N
   e^{i\mathbf{q}\cdot(\mathbf{r}_k(t)-\mathbf{r}_j(0))}\right\rangle

which itself is equal to the static structure factor when :math:`t=0`, i.e.

.. math::

   S(q) = F(q,0)




All four quantities are calculated in histogram form by binning into *bin*
bins in Fourier space with spacing :math:`\Delta q`, where
:math:`\Delta q` is the maximum of the set
:math:`\{2\pi/L_x,2\pi/L_y,2\pi/L_z\}` (:math:`L_x` being the size of the
box domain in the :math:`x`-axis, etc).

For *structure/sincos*, *bin* is the total number of bins that will be
used, i.e. it will calculate quantities for multiple :math:`q` values.

For *structure/sincos/single*, *bin* represents the single bin (and so
single :math:`q` value) that quantities will be calculated at. We refer to
this as the *single* version of this compute.

The *single* is included in this package because sometimes one is
interested in only a single value of :math:`q` when e.g. examining the intermediate
scattering function. Additionally, only the *single* version of this compute is
compatible with e.g. :doc:`fix ave/correlate <fix_ave_correlate>` since it outputs
a vector (instead of an array).


Output info
"""""""""""

*structure/sincos* calculates a global array in which the number of rows is
*bins* and the number of columns is 4. *structure/sincos/single* calculates a global
vector with 4 columns (each with a single entry).

For both cases, the first column is the
(solid angle) average value of the vector :math:`\mathbf{q}` in the
binned spherical (or circular, for 2D simulations) shells, so the true :math:`q`
value of a bin is given by :math:`q=E[\mathbf{q}]`
where :math:`E[]` indicates the solid angle average.

The second column is the count of :math:`\mathbf{q}` values which land in the
binned spherical (or circular, for 2D simulation) shells, which we will
call :math:`N_{counts}(q)`. 

The third and fourth columns compute the solid-angle averaged quantities

.. math::

   \zeta_1(q,t) = \sqrt{N_{counts}(q)}E\bigg[\sum_{j=1}^N \cos(\mathbf{q}\cdot\mathbf{r}_j(t))\bigg]

   \zeta_2(q,t) = \sqrt{N_{counts}(q)}E\bigg[\sum_{j=1}^N \sin(\mathbf{q}\cdot\mathbf{r}_j(t))\bigg]



The factors of :math:`\sqrt{N_{counts}(q)}` are necessary to ensure that squaring these
quantities will yield the correct counting :math:`N_{counts}(q)` of 
These two quantities are useful because

.. math::

   \sum_{j=1}^N\sum_{k=1}^N
   e^{i\mathbf{q}\cdot(\mathbf{r}_k(t)-\mathbf{r}_j(0))}
   =\sum_{j=1}^N
   \bigg(\cos(\mathbf{q}\cdot\mathbf{r}_k(t))
   +i\sin(\mathbf{q}\cdot\mathbf{r}_k(t))\bigg)
   \sum_{k=1}^N\bigg(\cos(\mathbf{q}\cdot\mathbf{r}_j(0))
   -i\sin(\mathbf{q}\cdot\mathbf{r}_j(0))\bigg)
   
which means that one can directly calculate the intermediate scattering
function and static structure factor since e.g.

.. math::

   F(q,t) = \left\langle\frac{1}{N}\big(\zeta_1(q,t)+i\zeta_2(q,t)\big)
   \big(\zeta_1(q,0)-i\zeta_2(q,0)\big)\right\rangle




These values can be used
by any command that uses a global values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.


The first column of array values will be in inverse distance
:dod:`units <units>`. The remaining three columns are dimensionless.
The first two columns
will remain unchanged throughout the simulation unless the simulation
box size is changing dynamically.

Restrictions
""""""""""""

Currently, this compute does not support distinguishing between different
types of particles (see e.g. :doc:`compute rdf <compute_rdf>`).

This compute requires a for loop to iterate through Fourier space. Therefore,
it scales with *bin* to the power of simulation dimension. This means the
compute is not particularly efficient and should be used sparingly
throughout the simulation.

This compute is part of the STRUCTURE-FACTOR package.  It is only
enabled if LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

The compute_structure_factor command does not work for triclinic cells.


Related commands
""""""""""""""""


:doc:`compute intermediatescattering <compute_intermediatescattering>`
:doc:`compute intermediatescattering/single <compute_intermediatescattering_single>`
:doc:`compute xrd <compute_xrd>`

Default
"""""""


