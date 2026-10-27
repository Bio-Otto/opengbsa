Generalized Born models
=======================

A Generalized Born (GB) model gives the polar part of the solvation free
energy. Each atom has an effective Born radius that measures how deeply it is
buried in the solute; the GB energy then combines the charges, the distances
and the Born radii. The models differ in how they compute the Born radii.
OpenGBSA provides the five GB models of Amber, set with
``analysis_settings.gb_model``:

.. list-table::
   :header-rows: 1
   :widths: 12 14 16 58

   * - ``gb_model``
     - Amber ``igb``
     - Radii set
     - Reference
   * - ``HCT``
     - 1
     - ``mbondi``
     - Hawkins GD, Cramer CJ, Truhlar DG (1995) *Chem Phys Lett* 246:122-129
   * - ``OBC1``
     - 2
     - ``mbondi2``
     - Onufriev A, Bashford D, Case DA (2000) *J Phys Chem B* 104:3712-3720
   * - ``OBC2``
     - 5
     - ``mbondi2``
     - Onufriev A, Bashford D, Case DA (2004) *Proteins* 55:383-394
   * - ``GBn``
     - 7
     - ``mbondi`` (``bondi`` in Amber)
     - Mongan J, Simmerling C, McCammon JA, Case DA, Onufriev A (2007)
       *J Chem Theory Comput* 3:156-169
   * - ``GBn2``
     - 8
     - ``mbondi3``
     - Nguyen H, Roe DR, Simmerling C (2013) *J Chem Theory Comput*
       9:2020-2034

**HCT** sums pairwise descreening contributions and tends to underestimate the
Born radii of buried atoms. **OBC1** and **OBC2** rescale the HCT integral to
correct this, with two published parameter sets; OBC2 is the default of
OpenGBSA and the most widely used model in MM/GBSA. **GBn** adds a correction
for the small regions between atoms ("necks") that are too narrow for water,
and **GBn2** refines it with new parameters fitted for proteins. All five
models were compared with Amber MMPBSA.py on the same frames
(:doc:`../testing`).

Atomic radii
------------

Each model was parameterised with its own set of intrinsic atomic radii, listed
above. OpenGBSA chooses them as follows (``gb_radii_set: auto``):

* An Amber ``.prmtop`` keeps the radii it contains, as Amber MMPBSA.py does.
  OpenGBSA warns when they differ from the set the model expects.
* GROMACS, CHARMM/NAMD, OpenMM and PDB input carry no GB radii. OpenGBSA
  assigns the model's set from each atom's element and the atom it is bonded
  to, by the rules of Amber's ``tleap``.
* GBn and GBn2 are defined only for radii between 1 and 2 Å. Radii outside
  this range are replaced by the model's own set, with a warning.

``gb_radii_set`` forces one set (``mbondi``, ``mbondi2``, ``mbondi3`` or
``bondi``) for every input.

The Born-radius cutoff
----------------------

An atom's Born radius depends on all other atoms of the molecule. OpenGBSA
computes it from every atom pair by default. Amber includes only the atoms
within 25 Å (``rgbmax``). ``rgbmax`` sets such a cutoff in OpenGBSA, and
``amber_compatible: true`` sets 25 Å. The effect of the cutoff is explained in
:doc:`mmgbsa_theory`.

Salt
----

Mobile ions screen the electrostatic interactions. With
``salt_concentration`` :math:`I` (in M) the GB energy includes the
Debye-Hückel screening of Srinivasan et al. (1999), as in Amber:

.. math::

   \kappa = 0.73 \times 50.33355 \sqrt{\frac{I}{\varepsilon_\mathrm{solvent} T}}\ \text{Å}^{-1}

At 0.15 M, :math:`\varepsilon_\mathrm{solvent} = 78.5` and 300 K this gives
:math:`\kappa = 0.093` Å⁻¹. ``salt_concentration: 0`` switches screening off.

Dielectric constants
--------------------

The solute dielectric is 1.0 and the solvent dielectric 78.5 by default
(``solute_dielectric``, ``solvent_dielectric``). With ``amber_compatible``
the solvent dielectric of MMPBSA.py is used: 80.0 when it evaluates the energy
with ``sander`` (CHARMM and chamber topologies) and 78.3 when it uses NAB.
