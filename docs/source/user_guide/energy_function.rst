The energy function
===================

Every frame of an OpenGBSA run gives three energies, those of the complex, the
receptor and the ligand, and the binding enthalpy is their difference
(:doc:`../theory/mmgbsa_theory`). The options on this page decide how each of
the three energies is computed. They all belong to ``analysis_settings``.

.. code-block:: yaml

    analysis_settings:
      gb_model: OBC2              # Generalized Born model
      gb_radii_set: auto          # intrinsic radii of the GB model
      salt_concentration: 0.15    # mol/L, Debye-Hueckel screening
      temperature: 300.0          # K, screening and entropy
      solute_dielectric: 1.0      # inside the molecules
      solvent_dielectric: 78.5    # water
      sa_model: LCPO              # non-polar solvation term
      lcpo_typing: element        # how LCPO parameters are assigned
      rgbmax: null                # Born-radius cutoff in A; null = none
      amber_compatible: false     # Amber MMPBSA.py conventions

These are the defaults. The rest of this page explains each option: what it
does, how OpenGBSA applies it, how much it changes a result and when to change
it. The table at the end of the page shows the effect of each option on two
bundled systems.

How the three energies are built
--------------------------------

For each of the complex, the receptor and the ligand OpenGBSA builds one OpenMM
``System``. The molecular mechanics terms come from the topology (Native Mode)
or from a force field (Coordinate Mode, :doc:`input_formats`). OpenGBSA then
changes the system in four ways:

#. The van der Waals interactions are moved from the ``NonbondedForce`` into a
   ``CustomNonbondedForce``, so that van der Waals and electrostatic energies can
   be read separately.
#. A Generalized Born force for the polar solvation term is built, a
   ``CustomGBForce`` from OpenMM's ``openmm.app.internal.customgbforces``, with
   the model, radii, dielectric constants and salt screening set by the options
   below. A GB force that the topology reader may have created is replaced, so
   every model gets the same settings in every input format.
#. The electrostatic energy is divided by ``solute_dielectric``.
#. A force for the non-polar solvation term is added (``sa_model``).

No interaction is cut off: every atom pair is included, as in Amber MMPBSA.py.
Each term sits in its own force group and is evaluated separately for every
frame. These are the columns ``delta_vdw``, ``delta_elec``, ``delta_gb`` and
``delta_sa`` of the per-frame output (:doc:`outputs`). Bonded terms are
evaluated as well, but they cancel exactly because receptor and ligand take
their coordinates from the complex.

.. _opt-gb_model:

``gb_model``
------------

:Default: ``OBC2``
:Values: ``HCT``, ``OBC1``, ``OBC2``, ``GBn``, ``GBn2``
:Amber: ``igb`` = 1, 2, 5, 7, 8 in MMPBSA.py's ``&gb`` namelist

The Generalized Born model computes the polar part of the solvation free
energy from the charges, the distances and an effective Born radius for every
atom. The five models differ in how they compute the Born radii
(:doc:`../theory/gb_models`).

*How it is applied.* OpenGBSA builds the matching OpenMM force:
``GBSAHCTForce``, ``GBSAOBC1Force``, ``GBSAOBC2Force``, ``GBSAGBnForce`` or
``GBSAGBn2Force``. Each atom gets its charge, its intrinsic radius
(``gb_radii_set``) and its screening factor. HCT, OBC1 and OBC2 take the
screening factors of the topology. GBn and GBn2 use the per-element values they
were fitted with, and GBn2 also its per-element :math:`\alpha`, :math:`\beta`
and :math:`\gamma`, with separate values for nucleic-acid residues. A residue
counts as a nucleic acid by its name: ``A``, ``C``, ``G``, ``U``, ``RA``,
``RC``, ``RG``, ``RU``, ``DA``, ``DC``, ``DG`` and ``DT``, each also with Amber's
suffix ``5``, ``3`` or ``N`` for terminal residues as in ``sander``, and the
CHARMM names ``ADE``, ``CYT``, ``GUA``, ``THY`` and ``URA``. On a
DNA-peptide test system the GB energy of every model agrees with ``sander``
within 0.22 kcal/mol out of 1,100.

*Effect.* The model changes the polar term of every complex, so it changes the
binding enthalpy more than any other option except ``solute_dielectric``. The
measured effects are in the table at the end of this page.

.. TODO numbers from option_effects.csv

*When to change it.* OBC2 is the most widely used model in MM/GBSA and the
default. GBn2 was fitted for proteins and nucleic acids and is a good choice for
protein-DNA and protein-RNA complexes. To compare with published Amber results,
use the model of that study. Compare results only within one model.

.. code-block:: yaml

    analysis_settings:
      gb_model: GBn2

.. _opt-gb_radii_set:

``gb_radii_set``
----------------

:Default: ``auto``
:Values: ``auto``, ``mbondi``, ``mbondi2``, ``mbondi3``, ``bondi``
:Amber: the radii in the ``.prmtop``, set in ``tleap`` with ``set default PBRadii``

Each GB model was parameterised with a set of intrinsic atomic radii: ``mbondi``
for HCT, ``mbondi2`` for OBC1 and OBC2, ``bondi`` for GBn and ``mbondi3`` for
GBn2.

*How it is applied.* With ``auto``, an Amber ``.prmtop`` keeps the radii it
contains, as Amber does. OpenGBSA warns when they are not the set of the
chosen model. GBn and GBn2 are defined only for radii between 1 and 2 Å; radii
outside that range are replaced by the model's set, with a warning. Input
without GB radii (GROMACS, CHARMM/NAMD, OpenMM and PDB files) receives the
model's set, assigned from each atom's element and bonding with the rules of
Amber's ``tleap``. Any other value assigns that set to every atom, whatever the
input. The decision is taken once, for the complex, and used for receptor and
ligand as well, so the three always have the same radii.

*Effect.* The quickstart topology already carries ``mbondi2`` radii. Other sets
change the result by the amounts in the table at the end of this page.

*When to change it.* Keep ``auto``. Set the radii explicitly only to use a
set that differs from the one in the topology, for example to follow a
published protocol.

.. _opt-salt_concentration:

``salt_concentration``
----------------------

:Default: ``0.15``
:Values: 0 to 2, in mol/L
:Amber: ``saltcon``

Mobile ions screen electrostatic interactions over distances longer than the
Debye length, about 8 Å at 0.15 M. OpenGBSA includes this screening in the GB
term with the Debye-Hückel approximation of Srinivasan et al. (1999), as Amber
does.

*How it is applied.* The ionic strength :math:`I` gives the screening
parameter

.. math::

   \kappa = 0.73 \times 50.33355 \sqrt{\frac{I}{\varepsilon_\mathrm{solvent}\,T}}\ \text{Å}^{-1}

with the solvent dielectric and the temperature of the run. The factor 0.73
accounts for the volume from which ions are excluded. :math:`\kappa` enters
the GB energy as :math:`e^{-\kappa f_\mathrm{GB}}/\varepsilon_\mathrm{solvent}`
(:doc:`../theory/mmgbsa_theory`). It does not change the gas-phase electrostatic
energy. ``salt_concentration: 0`` switches screening off. With
``amber_compatible: true`` OpenGBSA uses Amber's constant
:math:`\kappa = 0.73\sqrt{0.10806\,I}` Å⁻¹ instead, which, as in ``sander``,
depends on neither the temperature nor the solvent dielectric.

*Effect.* Screening mostly acts between atoms far apart. Receptor and ligand
atoms that interact are close, so salt changes the binding enthalpy of a small
ligand very little. It matters for highly charged complexes, such as protein
and DNA.

*When to change it.* Use the ionic strength of the experiment or of the
simulation. Use 0 to compare with results computed without salt.

.. _opt-temperature:

``temperature``
---------------

:Default: ``300.0``
:Values: larger than 0 and up to 1000, in K
:Amber: none; ``sander``'s salt screening does not depend on ``temp0``, and MMPBSA.py's normal-mode entropy is fixed at 298.15 K

The temperature enters the salt screening parameter :math:`\kappa` above and
the entropy estimates (:doc:`entropy`). No molecular mechanics term depends on
it.

*How it is applied.* :math:`\kappa` varies as :math:`T^{-1/2}`, so 310 K
instead of 300 K lowers it by 1.6 %. Without salt the temperature does not
change the binding enthalpy at all. With ``amber_compatible: true`` the
screening does not depend on the temperature, as in ``sander``. The entropy
estimates use it as :math:`k_\mathrm{B}T` and :math:`T\Delta S`, except the
normal-mode estimate with ``amber_compatible: true``, which always uses
MMPBSA.py's 298.15 K.

*When to change it.* Use the temperature of the simulation.

.. _opt-solute_dielectric:

``solute_dielectric``
---------------------

:Default: ``1.0``
:Values: a number larger than 0
:Amber: ``sander``'s ``intdiel``; MMPBSA.py fixes it at 1

The dielectric constant inside the molecules. It scales every electrostatic
interaction of the solute.

*How it is applied.* As Amber's ``intdiel``, it acts on two terms. The
gas-phase electrostatic energy, 1-4 pairs included, is divided by
:math:`\varepsilon_\mathrm{solute}`. The GB energy is multiplied by
:math:`1/\varepsilon_\mathrm{solute} - e^{-\kappa f_\mathrm{GB}}/\varepsilon_\mathrm{solvent}`.
OpenGBSA divides the charges of the ``NonbondedForce`` by
:math:`\sqrt{\varepsilon_\mathrm{solute}}` after the GB force has been built from
the unscaled charges.

*Effect.* A larger value weakens all electrostatic contributions. Because the
electrostatic interaction and the polar solvation largely cancel, the binding
enthalpy changes by the difference of two large terms, which can be large.

*When to change it.* Keep 1, the usual value in MM/GBSA and the value
MMPBSA.py uses. Values of 2 or 4 are sometimes used for binding sites with many
charged residues. Test such a choice against experimental data before relying
on it.

.. code-block:: yaml

    analysis_settings:
      solute_dielectric: 2.0

.. _opt-solvent_dielectric:

``solvent_dielectric``
----------------------

:Default: ``78.5``, or MMPBSA.py's value with ``amber_compatible``
:Values: a number larger than 0
:Amber: ``sander``'s ``extdiel``; MMPBSA.py fixes it at 80.0 or 78.3

The dielectric constant of water. It enters the GB energy and the salt
screening parameter.

*How it is applied.* The GB energy is multiplied by
:math:`1/\varepsilon_\mathrm{solute} - e^{-\kappa f_\mathrm{GB}}/\varepsilon_\mathrm{solvent}`.
With ``amber_compatible: true`` and no explicit value, OpenGBSA uses the value
that MMPBSA.py uses for the topology: 80.0 when MMPBSA.py evaluates the energy
with ``sander``, which it does for CHARMM and chamber topologies, and 78.3 when
it uses its NAB program, which it does otherwise. OpenGBSA counts a topology as
CHARMM-like when it is a ``.psf``, a chamber ``.prmtop`` or contains CMAP or
Urey-Bradley terms.

*Effect.* Small. Between 78.5 and 80 the factor :math:`1 - 1/\varepsilon_\mathrm{solvent}`
changes by 0.02 %.

*When to change it.* Rarely; keep the default.

.. _opt-sa_model:

``sa_model``
------------

:Default: ``LCPO``
:Values: ``LCPO``, ``ACE``
:Amber: ``surften`` (0.0072) and ``surfoff`` (0) of MMPBSA.py's GB calculation

The non-polar part of the solvation free energy, proportional to the
solvent-accessible surface area: :math:`G_\mathrm{nonpolar} = \gamma\,
\mathrm{SASA}`.

*How it is applied.* ``LCPO`` adds OpenMM's ``LCPOForce``: the surface area
from the linear combination of pairwise overlaps of Weiser et al. (1999) with a
probe radius of 1.4 Å, :math:`\gamma` = 0.0072 kcal/mol/Å² and no offset, as in
Amber's GB calculations. ``ACE`` adds OpenMM's analytical continuum
electrostatics approximation of the surface term, with :math:`\gamma` = 0.0054
kcal/mol/Å² and the element radii of OpenMM.

*Effect.* The two models give different surface areas, and ``ACE`` uses a
different :math:`\gamma`. See the table at the end of this page.

*When to change it.* Keep ``LCPO``, which is Amber's choice and the one used to
validate OpenGBSA. ``ACE`` is there for comparison with OpenMM's own implicit
solvent.

.. _opt-lcpo_typing:

``lcpo_typing``
---------------

:Default: ``element``, or ``amber`` with ``amber_compatible``
:Values: ``element``, ``amber``

LCPO has different parameters for different kinds of atoms: an sp3 carbon with
one heavy-atom neighbour, a carbonyl oxygen and so on.

*How it is applied.* ``element`` assigns them from each atom's element and its
bonded neighbours, which does not depend on how a force field names its atom
types. ``amber`` assigns them as ``sander`` does, from the element and the first
letters of the Amber atom-type name. The two agree for Amber force fields. They
differ for atom types that do not follow Amber's names, such as CHARMM types or
some ligand types.

*Effect.* Only on the non-polar term, and only where the two rules differ.

*When to change it.* Keep the default. ``amber_compatible`` selects ``amber``
so that the non-polar term matches MMPBSA.py exactly.

.. _opt-rgbmax:

``rgbmax``
----------

:Default: ``null`` (no cutoff), or 25 with ``amber_compatible``
:Values: ``null`` or a distance larger than 0, in Å
:Amber: ``rgbmax``

An atom's Born radius depends on how deeply it is buried, which is computed
from all other atoms. Amber includes only atoms within ``rgbmax``, 25 Å in
MMPBSA.py.

*How it is applied.* With a value, an atom :math:`j` contributes to the Born
radius of atom :math:`i` only if it lies within ``rgbmax`` of it, and an atom
that straddles the cutoff contributes only the part inside it, exactly as in
``sander``. The GB energy itself still includes every atom pair.

*Effect.* Pairs within the receptor or within the ligand are affected in the
same way in the complex and in the free molecules, and their effect largely
cancels. Pairs between receptor and ligand exist only in the complex, so the
cutoff does not cancel for them. The difference grows with the size of the
complex, from about 0.1 kcal/mol for small ligands to 16 kcal/mol for an
antibody-gp120 complex (:doc:`../theory/mmgbsa_theory`).

*When to change it.* Keep ``null``. The cutoff has no physical basis; it exists
to make Amber's calculation faster. Use ``amber_compatible`` to compare with
MMPBSA.py.

.. _opt-amber_compatible:

``amber_compatible``
--------------------

:Default: ``false``
:Values: ``true``, ``false``

Applies Amber MMPBSA.py's conventions in one switch, to compare OpenGBSA with
MMPBSA.py results on the same frames:

* ``rgbmax`` 25 Å;
* MMPBSA.py's solvent dielectric for the topology, 80.0 or 78.3;
* ``lcpo_typing: amber``;
* Amber's salt screening parameter, which does not depend on the temperature;
* for the normal-mode entropy, MMPBSA.py's energy function and temperature
  (:doc:`entropy`).

An option set explicitly in the configuration keeps its value. ``gb_model``,
``salt_concentration``, ``gb_radii_set`` and ``solute_dielectric`` are not
changed, so set them to the values of the MMPBSA.py input. With
``amber_compatible`` OpenGBSA reproduced MMPBSA.py to within 0.05 kcal/mol for
every validation system (:doc:`../testing`).

.. code-block:: yaml

    analysis_settings:
      gb_model: OBC2            # MMPBSA.py: igb = 5
      salt_concentration: 0.15  # MMPBSA.py: saltcon = 0.15
      amber_compatible: true

Options that were removed
-------------------------

``nonbonded_cutoff`` was removed in version 0.0.23. It applied a cutoff to
some interactions and not to others, and on GPUs it stopped the run. MM/GBSA
evaluates every atom pair, as MMPBSA.py does. Use ``rgbmax`` to limit the
Born-radius calculation.

Effect of each option
---------------------

.. TODO csv-table of docs/data/option_effects.csv after the code changes are merged

The table shows the mean binding enthalpy of the quickstart system
(:doc:`../quickstart`, 50 frames) when one option is changed from the
quickstart settings, and the mean of each term. It is produced by
``docs/scripts/option_effects.py``.
