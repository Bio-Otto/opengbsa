"""
LCPO parameter assignment as Amber MMPBSA.py does it (``lcpo_typing: amber``).

LCPO (Weiser, Shenkin and Still, J. Comput. Chem. 1999, 20, 217) gives every
atom a sphere radius and four coefficients that depend on its chemistry. OpenGBSA
assigns them by default from each atom's element and bonding (OpenMM's
``getLCPOParamsTopology``). Amber assigns them differently: ``sander``'s
``gbsa=1`` code, which MMPBSA.py uses for its non-polar term, chooses by element
and then by the *name* of the atom type (``O``, ``O2``, ``N3``, ``SH``, ``MG``,
``F``, ``Z*``), and gives any other element the sp2-carbon parameters. The rule
below follows ``Action_Surf::SetAtomLCPO`` in cpptraj, which states that it is
adapted from the ``gbsa=1`` method in sander's ``mdread.F90``.

For Amber force-field types the two assignments agree. For other types they do
not: a CHARMM carboxylate oxygen (``OC``) receives hydroxyl parameters, a CHARMM
lysine nitrogen (``NH3``) sp2-nitrogen parameters, and a GAFF chlorine (``cl``)
sp2-carbon parameters. The non-polar term then differs from Amber's by up to a
few tenths of a kcal/mol. ``amber_compatible`` applies this rule so that the
non-polar term matches MMPBSA.py on the same topology.
"""
from __future__ import annotations

import logging

try:
    from openmm.app.internal.lcpo import LCPO_PARAMETERS
except ImportError:  # pragma: no cover - OpenMM without LCPO support
    LCPO_PARAMETERS = None

log = logging.getLogger(__name__)


def amber_lcpo_key(element: int, type_name: str, n_heavy: int, n_total: int) -> str:
    """Key into ``LCPO_PARAMETERS`` chosen by Amber's ``gbsa=1`` rule.

    element: atomic number (0 if unknown); type_name: the atom-type name in the
    topology; n_heavy: bonded non-hydrogen neighbours; n_total: all bonded neighbours.
    """
    t = (type_name or "").upper()
    t0, t1 = t[:1], t[1:2]
    if element == 6:
        if n_total == 4:
            return {1: "C_sp3_1", 2: "C_sp3_2", 3: "C_sp3_3", 4: "C_sp3_4"}.get(n_heavy, "C_sp3_1")
        return {2: "C_sp2_2", 3: "C_sp2_3"}.get(n_heavy, "C_sp3_1")
    if element == 8:
        if t == "O":
            return "O_sp2_1"
        if t0 == "O" and t1 == "2":
            return "O_carboxylate"
        return {1: "O_sp3_1", 2: "O_sp3_2"}.get(n_heavy, "O_sp3_1")
    if element == 7:
        if t0 == "N" and t1 == "3":
            return {1: "N_sp3_1", 2: "N_sp3_2", 3: "N_sp3_3"}.get(n_heavy, "N_sp3_1")
        return {1: "N_sp2_1", 2: "N_sp2_2", 3: "N_sp2_3"}.get(n_heavy, "N_sp3_1")
    if element == 16:
        return "S_1" if (t0 == "S" and t1 == "H") else "S_2"
    if element == 15:
        return {3: "P_3", 4: "P_4"}.get(n_heavy, "P_3")
    if element == 1:
        return "none"
    if t0 == "Z":
        return "none"
    if t0 == "M" and t1 == "G":
        return "Mg"
    if t0 == "F":
        return "F"
    return "C_sp2_2"          # Amber: "carbon SA parms for unknown atom"


def amber_lcpo_params(topology, type_names):
    """Per-atom LCPO parameters (radius in Angstrom, P1, P2, P3, P4 in 1/Angstrom^2)
    for an OpenMM ``Topology``, assigned with Amber's rule from ``type_names``."""
    if LCPO_PARAMETERS is None:
        raise ImportError("OpenMM LCPO parameters are not available")
    atoms = list(topology.atoms())
    if len(type_names) != len(atoms):
        raise ValueError(f"{len(type_names)} atom-type names for {len(atoms)} atoms")
    z = [0 if a.element is None else a.element.atomic_number for a in atoms]
    n_heavy, n_total = [0] * len(atoms), [0] * len(atoms)
    for a1, a2 in topology.bonds():
        i, j = a1.index, a2.index
        n_total[i] += 1
        n_total[j] += 1
        if z[j] > 1:
            n_heavy[i] += 1
        if z[i] > 1:
            n_heavy[j] += 1
    return [LCPO_PARAMETERS[amber_lcpo_key(z[k], type_names[k], n_heavy[k], n_total[k])]
            for k in range(len(atoms))]


def type_names_from_structure(structure):
    """Atom-type names of a ParmEd ``Structure`` (as written in the topology file)."""
    try:
        names = [a.type for a in structure.atoms]
    except AttributeError:
        return None
    return names if all(isinstance(n, str) and n for n in names) else None


def type_names_from_forcefield(forcefield, topology):
    """Amber atom-type class of every atom, as an OpenMM ``ForceField`` assigns it.

    These are the names ``tleap`` writes into a topology built with the same force
    field. Uses OpenMM's template matching; returns None if that is unavailable."""
    if forcefield is None:
        return None
    try:
        from openmm.app.forcefield import ForceField
        data = ForceField._SystemData(topology)
        forcefield._matchAllResiduesToTemplates(data, topology, {}, False)
        return [forcefield._atomTypes[data.atomType[a]].atomClass for a in topology.atoms()]
    except Exception as exc:  # internal API; fall back to the element-based rule
        log.warning(f"Could not read atom-type classes from the force field ({exc}); "
                    f"Amber LCPO typing is not available for this input.")
        return None
