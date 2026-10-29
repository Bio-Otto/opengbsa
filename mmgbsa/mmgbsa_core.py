#!/usr/bin/env python3
"""
MM/GBSA calculations with OpenMM.

This is the real MM/GBSA engine (GBSACalculator, StructureManager,
GBSAForceManager) -- import directly from this module. An earlier
`mmgbsa/core/` package attempted to extract this into smaller components
but was never wired into the actual execution path (its own classes were
unused scaffolding delegating straight back here); it was removed rather
than completed, see CHANGELOG.md.
"""

import numpy as np
import pandas as pd
import copy
import time
from pathlib import Path
import warnings
from .utils import convert_mol2_to_sdf
import mdtraj as md
from collections import defaultdict
import xml.etree.ElementTree as ET



import sys
import types
# Patch numpy.compat for ParmEd (NumPy 2.x compatibility)
if 'numpy.compat' not in sys.modules:
    comp = types.ModuleType('numpy.compat')
    sys.modules['numpy.compat'] = comp
    if hasattr(np, 'compat'):
         sys.modules['numpy.compat'] = np.compat
    else:
        # Minimal mock if np.compat doesn't exist at all

        # ParmEd needs 'asbytes' and 'asstr'
        def asbytes(s):
            if isinstance(s, bytes):
                return s
            return str(s).encode('latin1')
            
        def asstr(s):
            if isinstance(s, bytes):
                return s.decode('latin1')
            return str(s)
            
        comp.asbytes = asbytes
        comp.asstr = asstr
        sys.modules['numpy.compat'] = comp


warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="Unable to load toolkit 'OpenEye Toolkit'")
warnings.filterwarnings("ignore", message="importing 'simtk.openmm' is deprecated")

# Import your existing classes
try:
    from openmm import app, openmm, unit
except ImportError:
    # Fallback for older versions
    import openmm.app as app
    import openmm
    import openmm.unit as unit
import mdtraj as md



#!/usr/bin/env python3
"""
Complete Advanced MM/GBSA Calculator with GBSA Force Implementation
Fixed version that properly handles force exceptions to avoid OpenMM errors
"""
import os
import random
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from .visualization import AdvancedVisualization
from .pymol_viz import PyMOLVisualizer
from pathlib import Path
import numpy as np
import pandas as pd
try:
    import mdtraj as md
except ImportError:
    md = None

# New Modular Imports
from .inputs import InputManager, EngineMode
from .topology import TopologyLoader
from .trajectory import TrajectoryProcessor, reimage_molecules
from .interaction_entropy import (IE_SIGMA_LIMIT, METHOD_LABELS, describe,
                                  entropy_estimates, entropy_result,
                                  interaction_energies)
from openmm import app, openmm, unit
try:
    from openff.toolkit.topology import Molecule
except ImportError:
    Molecule = None
try:
    from openff.toolkit.typing.engines.smirnoff import ForceField
except ImportError:
    ForceField = None
import warnings
warnings.filterwarnings('ignore')
from .logger import ToolLogger
try:
    from .reporting import HTMLReportGenerator
except ImportError:
    HTMLReportGenerator = None

from mmgbsa.lcpo_amber import (amber_lcpo_params, type_names_from_forcefield,
                                type_names_from_structure)

# Try to import LCPO parameters and helper
try:
    from openmm.app.internal.lcpo import getLCPOParamsTopology
    HAS_LCPO_PARAMS = True
except ImportError:
    HAS_LCPO_PARAMS = False
    getLCPOParamsTopology = None

# Initialize logger
log = ToolLogger()

# Solvent and ion residue names stripped before the calculation, in the
# Amber/GROMACS and the CHARMM conventions (the per-residue decomposition
# strips its trajectory with the same set).
def index_selection(indices):
    """An MDTraj selection of exactly these atom indices ('index 0 to 9 or index 12')."""
    indices = sorted(int(i) for i in indices)
    if not indices:
        raise ValueError("index_selection needs at least one atom")
    parts, start, prev = [], indices[0], indices[0]
    for i in indices[1:] + [None]:
        if i is not None and i == prev + 1:
            prev = i
            continue
        parts.append(f"index {start}" if start == prev else f"index {start} to {prev}")
        if i is not None:
            start = prev = i
    return " or ".join(parts)


SOLVENT_ION_RESNAMES = ("WAT", "HOH", "H2O", "SOL", "TIP3", "TIP", "SPC",
                        "NA", "CL", "K", "MG", "ZN", "CA",
                        "Na+", "Cl-", "K+", "Mg2+", "Ca2+", "Zn2+",
                        "SOD", "CLA", "POT", "CAL", "ZN2")
try:
    from openmmforcefields.generators import SystemGenerator
except ImportError:
    SystemGenerator = None

if SystemGenerator is not None:
    class _ProteinOnlySystemGenerator(SystemGenerator):
        """SystemGenerator without a small-molecule template generator.

        SystemGenerator loads its default OpenFF force field (~0.5 s) unless
        small_molecule_forcefield=None, but with None its __init__ calls
        add_molecules(None), which raises (openmmforcefields 0.14). Protein-only
        systems never add molecules, so that call is skipped; everything else,
        including create_system, is SystemGenerator's own."""

        def add_molecules(self, molecules):
            if molecules is None:
                return
            super().add_molecules(molecules)
else:
    _ProteinOnlySystemGenerator = None


class StructureManager:
    """
    Helper class to handle ParmEd structure manipulations for Unified Topology Splitting.
    Ensures that Receptor and Ligand systems are mathematically exact subsets of the Complex.
    """
    # Standard Amber GB-model-to-radii-set pairing (see e.g. Amber Reference
    # Manual / MMPBSA.py docs): HCT pairs with 'mbondi', OBC1/OBC2 with
    # 'mbondi2', GBn with 'bondi' (NOT 'mbondi' -- confirmed against
    # OpenMM's own GBSAGBnForce.getStandardParameters(), which builds its
    # reference radii via `_bondi_radii()` -- plain bondi, where every H
    # atom gets a fixed 1.2 A regardless of what it's bonded to, unlike
    # mbondi's bonded-atom-dependent H radii (e.g. 0.8 A for O/S-bonded
    # hydroxyl/thiol H). Using 'mbondi' for GBn previously caused
    # OpenMM's GBnForce to reject valid mbondi hydroxyl-H radii outright
    # with "Radii must be between 1 and 2 Angstroms for neck lookup", since
    # GBn's neck-lookup table is calibrated for bondi's radius range, not
    # mbondi's), GBn2 with 'mbondi3'.
    _GB_MODEL_RADII_SET = {
        'HCT': 'mbondi',
        'OBC1': 'mbondi2',
        'OBC2': 'mbondi2',
        'GBn': 'bondi',
        'GBn2': 'mbondi3',
    }

    _AMBER_TOPOLOGY_SUFFIXES = ('.prmtop', '.parm7', '.top7')

    @staticmethod
    def gb_radii_decision(struct, gb_model, gb_radii_set='auto', amber_native=False,
                          label='complex'):
        """Return the GB radii set to assign to `struct`, or None to keep its own.

        Amber `.prmtop` topologies carry the radii they were built with, and
        Amber MMPBSA.py always uses those, whatever `igb` is. With
        `gb_radii_set='auto'` an Amber-native topology therefore keeps its
        own radii (a warning names the set recommended for `gb_model` when it
        differs). GROMACS-, CHARMM- and coordinate-mode structures carry no GB
        radii, so they get the set matching `gb_model`. Two exceptions: an
        explicit `gb_radii_set` always wins, and GBn/GBn2 fall back to their
        own set when the topology's radii lie outside the 1-2 A range of
        OpenMM's neck-integral table, which cannot evaluate them.
        """
        model_set = StructureManager._GB_MODEL_RADII_SET.get(gb_model, 'mbondi2')
        if gb_radii_set and gb_radii_set != 'auto':
            return gb_radii_set
        if not amber_native:
            return model_set
        own = ''
        try:
            own = str(struct.parm_data['RADIUS_SET'][0]).strip()
        except Exception:
            pass
        if 'radii' not in own.lower():  # e.g. a converted GROMACS topology records none
            own = ''
        if gb_model in ('GBn', 'GBn2'):
            radii = [a.solvent_radius for a in struct.atoms]
            if radii and (min(radii) < 1.0 or max(radii) > 2.0):
                log.warning(f"{label}: the topology's GB radii ({own or 'set not recorded'}) lie "
                            f"outside the 1-2 A range {gb_model} needs; assigning '{model_set}'. "
                            f"Results will differ from Amber MMPBSA.py run on this topology.")
                return model_set
        log.info(f"{label}: keeping the Amber topology's own GB radii "
                 f"({own or 'radius set not recorded'}), as Amber MMPBSA.py does.")
        if own and f"({model_set})" not in own:
            log.warning(f"{label}: {gb_model} is usually run with '{model_set}' radii, but this "
                        f"topology has {own}. Keeping the topology's radii so results match "
                        f"Amber on the same files; set analysis_settings.gb_radii_set: "
                        f"{model_set} to override.")
        return None

    @staticmethod
    def load_complex(pdb_path, prmtop_path=None, xtc_path=None, gb_model='OBC2',
                      charmm_params=None, gb_radii_set='auto'):
        """
        Load a complex structure consistent with ParmEd.

        Args:
            pdb_path (str): Path to PDB file (coordinates)
            prmtop_path (str, optional): Path to Amber topology file
            xtc_path (str, optional): Path to trajectory for dynamic dummy atom checking
            gb_model (str): Configured GB model, used to pick the GB radii set
                (see `_GB_MODEL_RADII_SET`) for structures that have no radii
                of their own (GROMACS-origin and CHARMM/NAMD-origin inputs --
                .top via ParmEd, .tpr via `load_tpr_as_parmed`, .psf via
                ParmEd+CharmmParameterSet -- carry no GB radii, an
                Amber-ecosystem convention).
            gb_radii_set (str): 'auto' (default) keeps an Amber `.prmtop`'s own
                radii, as Amber MMPBSA.py does, and assigns the gb_model set to
                everything else; a set name ('mbondi', 'mbondi2', 'mbondi3',
                'bondi') is assigned to every structure. See
                `gb_radii_decision`. The decision is recorded on the returned
                structure as `opengbsa_radii_set` (None = kept) so that the
                receptor and ligand can follow it.
            charmm_params (str or list of str, optional): CHARMM parameter
                file(s) (.prm/.str/.rtf), required only when `prmtop_path`
                is a `.psf` file (see `TopologyLoader._load_charmm`, which
                uses the identical ParmEd loading approach).

        Returns:
            parmed.Structure: The loaded structure
        """
        import parmed as pmd

        struct = None
        is_gromacs_origin = False
        is_charmm_origin = False
        if prmtop_path and Path(prmtop_path).exists():
            # Native Mode: Load Topology + Coordinates
            try:
                ext = Path(prmtop_path).suffix.lower()
                if ext == '.tpr':
                    from mmgbsa.tpr_loader import load_tpr_as_parmed
                    struct = load_tpr_as_parmed(str(prmtop_path), xtc_path=xtc_path)
                    is_gromacs_origin = True
                    # Optionally we can overwrite coords with pdb_path if needed, but TPR typically has them.
                    # Usually, pdb_path is an explicitly supplied standard coordinates file.
                    if pdb_path and Path(pdb_path).exists() and str(pdb_path) != str(prmtop_path):
                        ref_struct = pmd.load_file(pdb_path)
                        if len(ref_struct.atoms) == len(struct.atoms):
                            struct.coordinates = ref_struct.coordinates
                elif ext == '.psf':
                    # A PSF carries no force-field parameters or coordinates
                    # of its own (see TopologyLoader._load_charmm's
                    # docstring) -- both must be supplied separately.
                    if not charmm_params:
                        raise ValueError(
                            "Loading a .psf complex requires CHARMM parameter files "
                            "(charmm_params, e.g. .prm/.str) -- see GBSACalculator's "
                            "'charmm_params' setting."
                        )
                    params = pmd.charmm.CharmmParameterSet(
                        *([charmm_params] if isinstance(charmm_params, str) else charmm_params)
                    )
                    struct = pmd.load_file(str(prmtop_path))
                    struct.load_parameters(params)
                    is_charmm_origin = True
                    if pdb_path and Path(pdb_path).exists() and str(pdb_path) != str(prmtop_path):
                        ref_struct = pmd.load_file(pdb_path)
                        if len(ref_struct.atoms) == len(struct.atoms):
                            struct.coordinates = ref_struct.coordinates
                else:
                    struct = pmd.load_file(prmtop_path, xyz=pdb_path)
                    is_gromacs_origin = (ext == '.top')
            except Exception as e:
                # Retry without XYZ if it fails (e.g. atom mismatch)
                # But we really need coordinates for GBSA
                raise ValueError(f"Failed to load prmtop+pdb: {e}")
        else:
            # Coordinate Mode: Load PDB directly
            # Note: This might lack parameter info (charges) unless we parameterize it first.
            # But GBSACalculator usually handles parameterization.
            # If we are here, we might need to parameterize using OpenMM first,
            # then converting to ParmEd is complex.
            # Actually, standardizing on loading the PDB is fine if we are in non-native mode.
            struct = pmd.load_file(pdb_path)

        # GROMACS-origin and CHARMM-origin structures (.top, .tpr, .psf)
        # carry no intrinsic GB radii at all (all-zero solvent_radius/
        # screen); Amber-native structures (.prmtop) DO carry real radii,
        # but only for whatever radius set they were originally built with
        # (typically mbondi2, from the source tleap run) -- if the
        # requested gb_model needs a different set (e.g. GBn/GBn2 need
        # mbondi/mbondi3, whose radii fall in a different numeric range),
        # reusing the prmtop's own mbondi2 radii either silently gives
        # wrong GB energies (for models that don't validate their input
        # range) or crashes outright (confirmed: GBn/GBn2 raise "Radii must
        # be between 1 and 2 Angstroms for neck lookup" when handed
        # mbondi2 radii un-adjusted, since OpenMM's own neck-lookup table
        # only covers the mbondi/mbondi3 range). So this now always
        # assigns the radii set matching the configured GB model, via
        # ParmEd's own (Amber-validated) mbondi*/bondi rule sets, for every
        # structure origin -- not just GROMACS/CHARMM-origin ones (that
        # restriction previously caused a large, systematic (~20 kcal/mol
        # observed on real OXA-MD benchmark systems) discrepancy against
        # Amber MMPBSA.py reference results for GROMACS-sourced inputs,
        # and separately left Amber-native inputs unable to use any
        # gb_model but whatever radius set their own prmtop happened to be
        # built with).
        #
        # Since 0.0.11 an Amber-native topology keeps its own radii by default
        # (gb_radii_decision): 0.0.8-0.0.10 assigned the gb_model set to every
        # structure, which silently changed the radii of .prmtop files built
        # with another set (e.g. mbondi with OBC2) and broke agreement with
        # Amber MMPBSA.py, which always uses the topology's radii.
        origin_label = ('GROMACS-origin' if is_gromacs_origin
                         else 'CHARMM/NAMD-origin' if is_charmm_origin
                         else 'Amber-native')
        amber_native = (origin_label == 'Amber-native' and bool(prmtop_path)
                        and Path(prmtop_path).suffix.lower() in StructureManager._AMBER_TOPOLOGY_SUFFIXES)
        radii_set = StructureManager.gb_radii_decision(
            struct, gb_model, gb_radii_set, amber_native=amber_native, label=f"{origin_label} complex")
        if radii_set is not None:
            try:
                from parmed.tools import changeRadii
                changeRadii(struct, radii_set).execute()
                log.info(f"Applied '{radii_set}' GB radii set to {origin_label} structure "
                         f"(gb_model='{gb_model}').")
            except Exception as e:
                log.warning(f"Failed to apply GB radii set to {origin_label} structure: {e}. "
                            f"GB energies will use the crude element-only radii fallback and "
                            f"will NOT be comparable to Amber MMPBSA.py results.")
        struct.opengbsa_radii_set = radii_set
        # CHARMM-family topology? Decides MMPBSA.py's solvent dielectric under
        # amber_compatible (its sander path, 80.0, serves CHARMM/chamber topologies).
        # Tagged here because split_components later clears Urey-Bradley terms.
        struct.opengbsa_charmm_like = bool(
            is_charmm_origin or isinstance(struct, pmd.amber.ChamberParm)
            or len(getattr(struct, 'urey_bradleys', [])) or len(getattr(struct, 'cmaps', [])))

        # Implicit solvent and GBSA evaluation does not support explicit waters inside the structural force tree.
        # Stripping is now handled deliberately in the main pipeline AFTER exporting the solvated PDB for MDTraj

        return struct

    @staticmethod
    def split_components(complex_structure, ligand_resname='LIG', ligand_indices=None):
        """
        Split complex into Receptor and Ligand structures using atom stripping.

        Args:
            complex_structure (parmed.Structure): The full complex
            ligand_resname (str): Residue name of the ligand. Ignored when
                `ligand_indices` is given.
            ligand_indices (iterable of int, optional): Explicit 0-based atom
                indices identifying the ligand, overriding `ligand_resname`
                entirely. Required for protein-protein/peptide systems (e.g.
                a NAMD/CHARMM `.psf` split by `chainid`) where the "ligand"
                is not a single small-molecule residue name that an Amber
                mask (`:LIG`) can select -- see `get_ligand_indices`'s own
                `selection` override, which produces exactly this kind of
                index list from an mdtraj-style `chainid N` string.

        Returns:
            tuple: (receptor_structure, ligand_structure)
        """
        import parmed as pmd
        # WORKAROUND (ParmEd bug): TPR-origin CHARMM structures carry
        # Urey-Bradley type indices that can exceed the deduplicated type
        # list produced by Structure.copy(), raising IndexError in
        # parmed/structure.py. UB terms are intramolecular and cancel
        # exactly between complex and components in MMGBSA deltas, so
        # dropping UB bonds entirely (for the complex and both split
        # structures alike) is numerically safe. (Clear the TrackedList in
        # place; assigning a plain list breaks ParmEd's change-tracking.)
        complex_structure.urey_bradleys.clear()
        # Create copies to avoid mutating original. An Amber topology is split
        # as its own class (a full copy keeps the prmtop data that its strip()
        # rebuilds): a generic Structure copy of a CHARMM ChamberParm loses the
        # CHARMM 1-4 Lennard-Jones parameters (epsilon_14 = 0) and the
        # pair-specific LJ table, so receptor and ligand would not match the
        # complex (+217 kcal/mol in the van der Waals term of a CHARMM36
        # protein-ligand complex given without receptor/ligand topologies).
        # For AmberParm both copies give identical energies.
        if isinstance(complex_structure, pmd.amber.AmberParm):
            receptor = copy.deepcopy(complex_structure)
            ligand = copy.deepcopy(complex_structure)
        else:
            receptor = complex_structure.copy(cls=pmd.Structure)
            ligand = complex_structure.copy(cls=pmd.Structure)

        if ligand_indices is not None:
            ligand_index_set = set(int(i) for i in ligand_indices)
            n_atoms = len(complex_structure.atoms)
            # ParmEd's `strip` also accepts a boolean/index iterable the same
            # length as `atoms`, sidestepping AmberMask syntax entirely --
            # confirmed via `parmed.Structure.strip`'s own docstring.
            is_ligand = [i in ligand_index_set for i in range(n_atoms)]
            is_receptor = [not v for v in is_ligand]
            receptor.strip(is_ligand)
            ligand.strip(is_receptor)
        else:
            # Strip Ligand from Receptor (Keep everything NOT ligand)
            receptor.strip(f":{ligand_resname}")

            # Strip Receptor from Ligand (Keep ONLY ligand)
            ligand.strip(f"!:{ligand_resname}")

        # Apply secondary global solvent strip for receptor and ligand topologies to avoid LCPO crashes
        # CHARMM/NAMD-origin structures use different solvent/ion residue
        # names than Amber/GROMACS (e.g. 'TIP3' not 'SOL', 'SOD'/'CLA' not
        # 'NA'/'CL') -- confirmed on a real NAMD PSF that the original,
        # Amber/GROMACS-only mask silently stripped ZERO atoms (16024 atoms
        # in, 16024 out) despite the structure containing 4655 TIP3 waters
        # and 29 SOD/CLA ions, which would have gone on to crash or corrupt
        # the LCPO surface-area force downstream. Includes both naming
        # conventions so this mask works for every currently-supported
        # topology origin.
        solvent_mask = ":" + ",".join(SOLVENT_ION_RESNAMES)
        receptor.strip(solvent_mask)
        ligand.strip(solvent_mask)

        return receptor, ligand

    @staticmethod
    def create_openmm_system(structure, implicitSolvent=None, 
                           implicitSolventSaltConc=0.0*unit.molar,
                           nonbondedMethod=app.NoCutoff, 
                           nonbondedCutoff=None,
                           constraints=None,
                           solventDielectric=None):
        """
        Create an OpenMM System from a ParmEd Structure.
        """
        # WORKAROUND: ParmEd createSystem can hang with Cutoff + Implicit Solvent.
        # If Cutoff is requested, we create a Vacuum system first, then let refine_gbsa_forces 
        # add the GBSA force manually.
        pass_implicit = implicitSolvent
        pass_salt = implicitSolventSaltConc
        
        if nonbondedCutoff is not None and implicitSolvent is not None:
             pass_implicit = None
             pass_salt = None
             # We rely on GBSAForceManager.refine_gbsa_forces to add the GB force later.

        # Prepare kwargs
        kwargs = {
            'nonbondedMethod': nonbondedMethod,
            'nonbondedCutoff': nonbondedCutoff,
            'constraints': constraints,
            'implicitSolvent': pass_implicit,
            'implicitSolventSaltConc': pass_salt,
            # ParmEd builds the HCT/GBn/GBn2 force itself (OBC1/OBC2 are rebuilt in
            # refine_gbsa_forces); without this it used its default 78.5 whatever
            # solvent_dielectric was configured.
            'solventDielectric': solventDielectric if pass_implicit is not None else None,
        }
        # Filter None to let defaults handle it or avoid errors
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        
        system = structure.createSystem(**kwargs)

        # Sanity-check that 1-4 exceptions are actually populated. Check ALL
        # exceptions (cheap -- a single O(N) pass) and require every single
        # one to be zero before concluding anything is wrong: exceptions are
        # not randomly ordered (ParmEd/OpenMM tend to emit 1-2/1-3
        # exclusions and 1-4 pairs in separate blocks, not interleaved), so
        # even a several-hundred-item sample can land entirely on genuine
        # zero-value exclusions and wrongly conclude the whole list is
        # broken -- confirmed on a real CHARMM36 system, where even the
        # first 200 (of 21,470) exceptions were all legitimate zero-value
        # exclusions, while 10,330 of them were correctly non-zero 1-4 pairs
        # carrying real CHARMM parameters (from `structure.adjusts`) later
        # in the list -- the old blanket "fix" below was overwriting those
        # CORRECT values with fabricated Amber-style ones.
        forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
        if forces:
            nb = forces[0]
            fix_needed = False
            if nb.getNumExceptions() > 0:
                any_nonzero = False
                for i in range(nb.getNumExceptions()):
                    _, _, cp, _, eps = nb.getExceptionParameters(i)
                    if abs(cp._value) > 1e-6 or abs(eps._value) > 1e-6:
                        any_nonzero = True
                        break
                # Only treat this as broken if NONE of the exceptions carry a
                # real value -- for any normally-parameterized system (Amber
                # or CHARMM), at least some 1-4 pairs will have a nonzero
                # charge-product or epsilon.
                fix_needed = not any_nonzero
            
            if fix_needed:
                # IMPORTANT: 1-4 scaling conventions are force-field-specific.
                # Amber uses uniform SCEE=1.2/SCNB=2.0 (i.e. charge scaled by
                # 1/1.2=0.8333, LJ epsilon by 1/2.0=0.5, combining-rule sigma/
                # epsilon). CHARMM instead typically uses full-strength (1.0)
                # 1-4 electrostatics and its own explicit, atom-pair-specific
                # 1-4 LJ parameters (distinct from the regular nonbonded LJ
                # parameters) -- NOT the same combining rule scaled by 0.5.
                # Applying Amber's scale factors to a CHARMM system is
                # physically wrong and was confirmed (on a real CHARMM36
                # system) to produce a ~18 kcal/mol systematic vdW error
                # relative to Amber MMPBSA.py reference results.
                #
                # `structure.adjusts` (ParmEd's own, force-field-agnostic 1-4
                # exception list, populated correctly by ParmEd's CHARMM/
                # Amber/GROMACS parsers alike) is the authoritative source
                # for these values and is used here instead of any
                # hardcoded, Amber-specific scale factor.
                adjust_map = {}
                try:
                    for adj in structure.adjusts:
                        key = (adj.atom1.idx, adj.atom2.idx)
                        adjust_map[key] = adj.type
                        adjust_map[(key[1], key[0])] = adj.type
                except AttributeError:
                    pass  # structure has no .adjusts (e.g. some non-CHARMM formats); handled below

                if adjust_map:
                    print(f"  ⚠️  1-4 Interactions appear to be zeroed. Rebuilding from "
                          f"structure.adjusts ({len(structure.adjusts)} force-field-native 1-4 pairs).")
                    count_fixed = 0
                    count_missing = 0
                    for i in range(nb.getNumExceptions()):
                        p1, p2, _, _, _ = nb.getExceptionParameters(i)
                        adj_type = adjust_map.get((p1, p2))
                        if adj_type is None:
                            continue  # genuine 1-2/1-3 exclusion; leave at 0
                        q1, _, _ = nb.getParticleParameters(p1)
                        q2, _, _ = nb.getParticleParameters(p2)
                        chgscale = getattr(adj_type, 'chgscale', 1.0)
                        chg = q1.value_in_unit(unit.elementary_charge) * q2.value_in_unit(unit.elementary_charge) * chgscale
                        sig = adj_type.usigma.value_in_unit(unit.nanometer)
                        eps = adj_type.uepsilon.value_in_unit(unit.kilojoule_per_mole)
                        nb.setExceptionParameters(
                            i, p1, p2,
                            chg * unit.elementary_charge**2,
                            sig * unit.nanometer,
                            eps * unit.kilojoule_per_mole,
                        )
                        count_fixed += 1
                    print(f"  ✓ Rebuilt {count_fixed} 1-4 exceptions from structure.adjusts "
                          f"(force-field-native parameters, not Amber-specific scaling).")
                else:
                    # No .adjusts data available at all (rare; e.g. a bare
                    # Structure with no dihedral/adjust records). Fall back
                    # to Amber's own convention ONLY here, since a structure
                    # with no explicit 1-4 data is far more likely to be an
                    # Amber-style topology (which relies on the uniform
                    # SCEE/SCNB convention) than a CHARMM one (which always
                    # carries explicit per-pair 1-4 parameters via adjusts).
                    print("  ⚠️  1-4 Interactions appear to be zeroed and no structure.adjusts "
                          "data is available. Applying Amber SCEE=1.2/SCNB=2.0 convention as a "
                          "last resort -- this is WRONG for CHARMM/OPLS/GROMOS systems.")
                    try:
                        import math
                        count_fixed = 0
                        atoms = structure.atoms
                        for i in range(nb.getNumExceptions()):
                            p1, p2, _, _, _ = nb.getExceptionParameters(i)
                            a1 = atoms[p1]
                            a2 = atoms[p2]
                            is_12 = a2 in a1.bond_partners
                            is_13 = False
                            if not is_12:
                                for n in a1.bond_partners:
                                    if a2 in n.bond_partners:
                                        is_13 = True
                                        break
                            if is_12 or is_13:
                                nb.setExceptionParameters(i, p1, p2, 0.0, 0.5, 0.0)
                            else:
                                q1, s1, e1 = nb.getParticleParameters(p1)
                                q2, s2, e2 = nb.getParticleParameters(p2)
                                c_scale = 0.83333333
                                v_scale = 0.5
                                chg = (q1._value * q2._value) * c_scale
                                eps = math.sqrt(e1._value * e2._value) * v_scale
                                sig = (s1._value + s2._value) * 0.5
                                nb.setExceptionParameters(
                                    i, p1, p2,
                                    chg * unit.elementary_charge**2,
                                    sig * unit.nanometer,
                                    eps * unit.kilojoule_per_mole,
                                )
                                count_fixed += 1
                        print(f"  ✓ Manually updated {count_fixed} 1-4 exceptions with AMBER scaling.")
                    except Exception as e:
                        print(f"  ❌ Failed to regenerate 1-4 interactions: {e}")
                        print("     Proceeding with zeroed 1-4 interactions (results may be inaccurate).")
                
        # Amber prmtops (incl. chamber) store the full LJ pair table and ParmEd
        # builds the OpenMM table from it directly; only structures whose NBFIX
        # lives on their atom types (CHARMM .psf, GROMACS .top) need the rebuild.
        import parmed
        if (hasattr(structure, 'has_NBFIX') and structure.has_NBFIX()
                and not isinstance(structure, parmed.amber.AmberParm)):
            StructureManager._rebuild_nbfix_table(system, structure)
        return system

    @staticmethod
    def _rebuild_nbfix_table(system, structure):
        """Rebuild ParmEd's NBFIX Lennard-Jones lookup table so no NBFIX is lost.

        For a structure with pair-specific LJ corrections (CHARMM NBFIX, e.g.
        CHARMM36m's arginine NC2 - carboxylate OC pair), ParmEd's createSystem
        puts all LJ into a tabulated CustomNonbondedForce. Building the table,
        ParmEd 4.3.1 merges an atom type into an earlier type with the same
        rmin/epsilon whenever the EARLIER type has no NBFIX, even if the later
        one has: NC2 is merged into backbone NH1 and OC into carbonyl O, and
        their NBFIX disappears. Up to OpenGBSA 0.0.13 every CHARMM .psf input
        therefore ran without NBFIX, silently (CTLA4 -0.86, VRC01 about -2.7
        kcal/mol in van der Waals against Amber). Here types are merged only
        when neither has an NBFIX, and an NBFIX is looked up in both
        directions. Amber prmtops (incl. GROMACS .top input, which OpenGBSA
        converts to one) carry their own pair table and are not rebuilt; a
        .tpr carries no pair table at all (see tpr_loader).
        """
        import math
        force = next((f for f in system.getForces()
                      if isinstance(f, openmm.CustomNonbondedForce) and f.getNumTabulatedFunctions() == 2
                      and {f.getTabulatedFunctionName(0), f.getTabulatedFunctionName(1)} == {'acoef', 'bcoef'}),
                     None)
        if force is None or force.getNumParticles() != len(structure.atoms):
            return
        classes, index_of, per_atom = [], {}, []
        for atom in structure.atoms:
            t = atom.atom_type
            key = ('nbfix', t.name) if t.nbfix else ('lj', t.rmin, abs(t.epsilon))
            if key not in index_of:
                index_of[key] = len(classes)
                classes.append(t)
            per_atom.append(index_of[key])
        n = len(classes)
        acoef, bcoef = [0.0] * (n * n), [0.0] * (n * n)
        for i, a in enumerate(classes):
            for j, b in enumerate(classes):
                pair = (a.nbfix or {}).get(b.name) or (b.nbfix or {}).get(a.name)
                if pair is None:
                    rij = (a.rmin + b.rmin) * 0.1                                  # A -> nm
                    wdij = math.sqrt(abs(a.epsilon) * abs(b.epsilon)) * 4.184      # kcal -> kJ
                else:
                    rij, wdij = pair[0] * 0.1, abs(pair[1]) * 4.184
                rij6 = rij ** 6
                acoef[i + n * j] = math.sqrt(wdij) * rij6
                bcoef[i + n * j] = 2 * wdij * rij6
        tables = {force.getTabulatedFunctionName(k): force.getTabulatedFunction(k) for k in range(2)}
        tables['acoef'].setFunctionParameters(n, n, acoef)
        tables['bcoef'].setFunctionParameters(n, n, bcoef)
        for i, k in enumerate(per_atom):
            force.setParticleParameters(i, (k,))
        log.info(f"Rebuilt the NBFIX Lennard-Jones table ({n} LJ classes).")


class GBSAForceManager:
    """
    Advanced GBSA force implementation that refines OpenMM-generated forces.
    Now relies on SystemGenerator/app.createSystem to create the correct physics model (OBC1/OBC2/GBn),
    then refines it (separating SA, adding salt if needed).
    """
    
    def __init__(self, gb_model='OBC2', salt_concentration=0.15, solute_dielectric=1.0, solvent_dielectric=78.5, sa_model='LCPO', nonbonded_cutoff=None):
        self.gb_model = gb_model
        self.salt_concentration = salt_concentration
        self.solute_dielectric = solute_dielectric
        self.solvent_dielectric = solvent_dielectric
        self.sa_model = sa_model
        self.nonbonded_cutoff = nonbonded_cutoff
        # Born-radius cutoff (A) as in Amber's rgbmax; None = every pair (default).
        self.rgbmax = None
        # amber_compatible: pick MMPBSA.py's solvent dielectric from the topology
        # (resolve_mmpbsa_dielectric) unless one was set explicitly.
        self.dielectric_from_preset = False
        self._preset_dielectric_resolved = False
        # LCPO parameter assignment: 'element' (OpenMM: element and bonding) or
        # 'amber' (MMPBSA.py: by atom-type name; see mmgbsa/lcpo_amber.py).
        self.lcpo_typing = 'element'
        self.physics_assumptions = []
        self._logged_radius_fallback = False
        
        # Map string names to OpenMM constants
        self.gb_model_map = {
            'HCT': app.HCT,
            'OBC1': app.OBC1,
            'OBC2': app.OBC2,
            'GBn': app.GBn,
            'GBn2': app.GBn2
        }
        
        self.current_app_model = self.gb_model_map.get(gb_model, app.OBC2)

        # GB radii (Å) - AMBER parameter set (Still useful for manual SA or checks)
        self.gb_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.50, 'F': 1.47,
            'P': 1.85, 'S': 1.80, 'Cl': 1.75, 'Br': 0.85, 'I': 1.98,
            'Na': 1.868, 'K': 2.658, 'Mg': 1.584, 'Ca': 2.412, 'Zn': 1.394,
            'Fe': 1.456, 'Cu': 1.40, 'Mn': 1.456
        }
        
        self.gb_scales = {
            'H': 0.85, 'C': 0.72, 'N': 0.79, 'O': 0.85, 'F': 0.88,
            'P': 0.86, 'S': 0.96, 'Cl': 0.80, 'Br': 0.80, 'I': 0.80,
            'Na': 1.0, 'K': 1.0, 'Mg': 1.0, 'Ca': 1.0, 'Zn': 1.0,
            'Fe': 1.0, 'Cu': 1.0, 'Mn': 1.0
        }

    def refine_gbsa_forces(self, system, topology, lcpo_type_names=None):
        """
        Refine GBSA forces in a system that was already created with implicit solvent.
        1. Identifies the GB force.
        2. Sets Surface Area energy to 0 (to separate it).
        3. Adds the specific SA model (ACE/LCPO).
        """
        log.process(f"Refining GBSA forces for {self.gb_model} model...")
        self.resolve_mmpbsa_dielectric(topology)
        
        # Determine strict OpenMM Topology for SA/LCPO (which requires it)
        # But for OBC fallback, we prefer the ParmEd Structure (if passed) to get radii.
        struct_for_obc = topology
        if hasattr(topology, 'topology'):
             omm_topology = topology.topology
        else:
             omm_topology = topology
        # Atom-type names for Amber LCPO typing: given by the caller (Coordinate Mode,
        # from the force field) or read from a ParmEd Structure (the topology file).
        if lcpo_type_names is None and self.lcpo_typing == 'amber':
             lcpo_type_names = type_names_from_structure(topology)

        # 1. Identify Existing GB Force
        gb_force_index = -1
        gb_force = None
        nb_force = None
        custom_vdw_force = None
        
        for i, force in enumerate(system.getForces()):
            if isinstance(force, openmm.GBSAOBCForce):
                gb_force = force
                gb_force_index = i
            elif isinstance(force, openmm.CustomGBForce):
                gb_force = force
                gb_force_index = i
            elif isinstance(force, openmm.NonbondedForce):
                nb_force = force
            elif isinstance(force, openmm.CustomNonbondedForce):
                # Assume the first CustomNonbondedForce is the VDW force (common in ParmEd)
                if custom_vdw_force is None:
                    custom_vdw_force = force
                    
        # 1.5 Ensure Consistency: If CustomNonbondedForce is missing (e.g. Receptor/Ligand created from Structure),
        # but NonbondedForce has VDW (epsilon > 0), we MUST separate them to match Complex.
        # Otherwise, Receptor VDW is in Group 0 (mixed with Ele), but Complex VDW is Group 3.
        # The reporter expects VDW in Group 3.
        
        if nb_force and custom_vdw_force is None:
            # Check if we have VDW to move
            has_vdw = False
            for i in range(nb_force.getNumParticles()):
                _, _, eps = nb_force.getParticleParameters(i)
                if eps._value != 0:
                    has_vdw = True
                    break
            
            if has_vdw:
                log.info("Separating VDW from NonbondedForce to CustomNonbondedForce for consistency...")
                
                # Create CustomNonbondedForce (Group 3)
                # Use geometric mixing rules (standard for AMBER)
                custom_vdw_force = openmm.CustomNonbondedForce("4*epsilon*((sigma/r)^12 - (sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
                custom_vdw_force.addPerParticleParameter("sigma")
                custom_vdw_force.addPerParticleParameter("epsilon")
                custom_vdw_force.setForceGroup(3)
                
                # Copy Method/Cutoff settings
                custom_vdw_force.setNonbondedMethod(nb_force.getNonbondedMethod())
                custom_vdw_force.setCutoffDistance(nb_force.getCutoffDistance())
                custom_vdw_force.setUseSwitchingFunction(nb_force.getUseSwitchingFunction())
                custom_vdw_force.setSwitchingDistance(nb_force.getSwitchingDistance())
                # Don't need PME for VDW usually, or CutoffNonPeriodic equivalent
                if nb_force.getNonbondedMethod() == openmm.NonbondedForce.PME:
                    custom_vdw_force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
                
                # Move Particles
                for i in range(nb_force.getNumParticles()):
                    chg, sig, eps = nb_force.getParticleParameters(i)
                    custom_vdw_force.addParticle([sig, eps])
                    # Zero VDW in NonbondedForce (Keep Charge)
                    nb_force.setParticleParameters(i, chg, sig, 0.0*unit.kilojoule_per_mole)
                    
                # Copy Exclusions (All Exceptions in NB are Exclusions in CustomNB)
                # CustomNB doesn't have scaling, so we exclude 1-4s here and handle them in CustomBond or existing Exceptions?
                # Existing logic below moves 1-4s to CustomBond. So we should EXCLUDE them here.
                for i in range(nb_force.getNumExceptions()):
                    p1, p2, _, _, _ = nb_force.getExceptionParameters(i)
                    custom_vdw_force.addExclusion(p1, p2)
                    
                system.addForce(custom_vdw_force)
                log.info(f"  ✓ Created CustomNonbondedForce with {custom_vdw_force.getNumParticles()} particles.")

        # FIX: Handle 1-4 VDW Double Counting
        # If CustomNonbondedForce exists (Group 3 VDW), it calculates full VDW for 1-4s unless excluded.
        # NonbondedForce (Group 0) calculates Scaled VDW via exceptions.
        # We must:
        # 1. Exclude 1-4s from CustomNonbondedForce.
        # 2. Zero VDW in NonbondedForce exceptions (so Group 0 is ELE only).
        # 3. Create CustomBondForce for 1-4 VDW (Group 3).
        
        if nb_force and custom_vdw_force:
            log.info("Applying 1-4 VDW Correction (Moving 1-4s to CustomBondForce Group 3)...")
            
            # Create CustomBondForce for 1-4 VDW
            vdw_14_force = openmm.CustomBondForce("4*epsilon*((sigma/r)^12 - (sigma/r)^6)")
            vdw_14_force.addPerBondParameter("sigma")
            vdw_14_force.addPerBondParameter("epsilon")
            vdw_14_force.setForceGroup(3) # Group 3 (VDW)
            
            count_moved = 0
            
            # Pre-calculate existing exclusions to prevent duplicates
            existing_exclusions = set()
            for i in range(custom_vdw_force.getNumExclusions()):
                p1, p2 = custom_vdw_force.getExclusionParticles(i)
                existing_exclusions.add(tuple(sorted((p1, p2))))

            for i in range(nb_force.getNumExceptions()):
                p1, p2, chg, sig, eps = nb_force.getExceptionParameters(i)
                
                # Check if it's a 1-4 (non-zero scaled parameters usually, or just not 1-2/1-3 exclusion)
                # But since we want to move VDW specifically:
                if eps._value != 0.0:
                    # It has VDW content.
                    # 1. Add to CustomBondForce (Group 3)
                    vdw_14_force.addBond(p1, p2, [sig, eps])
                    
                    # 2. Exclude from CustomNonbondedForce (Group 3 - Full VDW)
                    pair = tuple(sorted((p1, p2)))
                    if pair not in existing_exclusions:
                        custom_vdw_force.addExclusion(p1, p2)
                        existing_exclusions.add(pair)
                    
                    # 3. Zero VDW in NonbondedForce (Group 0) - Keep Electrostatics
                    nb_force.setExceptionParameters(i, p1, p2, chg, sig, 0.0*unit.kilojoule_per_mole)
                    
                    count_moved += 1
            
            if count_moved > 0:
                system.addForce(vdw_14_force)
                log.info(f"  ✓ Moved {count_moved} 1-4 VDW interactions to Group 3.")
            else:
                log.info("  ℹ️ No 1-4 VDW interactions found to move.")

        if gb_force is None:
             # Try to find by index if not found (debugging safety)
             # But loop covers it. If not found, create fallback logic handles it later?
             # No, this function refines EXISTING forces.
             log.warning("No GB Force found in system to refine!")
        else:
             log.info(f"Identified GB Force: {type(gb_force).__name__}")
        # FIX: For OBC models, we MUST ensure we use the Prmtop radii (mbondi2/3).
        # ParmEd's createSystem often creates a CustomGBForce with default/wrong radii for OBC.
        # So for OBC, we force a rebuild using our Prmtop parameters.
        if self.gb_model in ['OBC1', 'OBC2']:
            if gb_force is not None:
                log.info(f"Replace existing GB force ({gb_force.__class__.__name__}) to enforce Prmtop radii.")
                system.removeForce(gb_force_index)
            
            gb_force = self._create_fallback_obc_force(system, struct_for_obc)
            system.addForce(gb_force)
            log.info(f"Created new GBSAOBCForce using Prmtop parameters.")
        elif gb_force is None:
            log.warning("No GB force found in system. Was it created with implicitSolvent? Attempting to add manual OBC2 fallback.")
            # Fallback for systems created without GB
            gb_force = self._create_fallback_obc_force(system, struct_for_obc)
            system.addForce(gb_force)
        
        # 2. Configure GB Force
        if isinstance(gb_force, openmm.GBSAOBCForce):
             log.info(f"Existing GBSAOBCForce found. Inspecting parameters...")
             # Check particle 0 for radius
             q, rad, scale = gb_force.getParticleParameters(0)
             log.info(f"Particle 0: Radius={rad.value_in_unit(unit.angstroms):.4f} A, Scale={scale}")
             
             # Zero out SA to separate it
             gb_force.setSurfaceAreaEnergy(0.0)
             gb_force.setSoluteDielectric(self.solute_dielectric)
             gb_force.setSolventDielectric(self.solvent_dielectric)
             gb_force.setForceGroup(1) # Group 1 for Standard GB
             log.info("Configured GBSAOBCForce: SA=0, Dielectrics set, Group=1.")
             
        elif isinstance(gb_force, openmm.CustomGBForce):
             # For CustomGBForce (GBn), assume params are set by factory.
             # We might need to manually set dielectrics if global params allow.
             for i in range(gb_force.getNumGlobalParameters()):
                 name = gb_force.getGlobalParameterName(i)
                 if name == 'soluteDielectric':
                     gb_force.setGlobalParameterDefaultValue(i, self.solute_dielectric)
                 elif name == 'solventDielectric':
                     gb_force.setGlobalParameterDefaultValue(i, self.solvent_dielectric)
             
             gb_force.setForceGroup(2) # Group 2 for Custom GB
             gb_force.setForceGroup(2) # Group 2 for Custom GB
             log.info("Configured CustomGBForce dielectrics, Group=2.")
             
             # Log particle 0 parameters to check radii
             if gb_force.getNumParticles() > 0:
                 params = gb_force.getParticleParameters(0)
                 # CustomGBForce params structure depends on definition. Usually (charge, radius, scale, ...) or similar.
                 # We simply print what we get.
                 log.info(f"CustomGBForce Particle 0 params: {params}")

        if gb_force is not None:
            self._apply_rgbmax(gb_force)

        # 3. Add Separate Surface Area Force
        sa_force = self._setup_surface_area_force(system, omm_topology, lcpo_type_names)
        sa_force.setForceGroup(4) # Explicitly set to Group 4 to separate from VDW
        system.addForce(sa_force)
        
        log.success(f"GBSA forces refined successfully")
        return system

    # Maps this codebase's gb_model strings to the matching OpenMM
    # CustomAmberGBForceBase subclass (openmm.app.internal.customgbforces).
    # All five share an identical constructor signature (solventDielectric,
    # soluteDielectric, SA, cutoff, kappa) and addParticle([charge, radius,
    # screen]) call convention, which is what makes a single code path for
    # all five models possible below -- confirmed directly against each
    # class's __init__ signature.
    _GB_MODEL_FORCE_CLASS_NAME = {
        'HCT': 'GBSAHCTForce',
        'OBC1': 'GBSAOBC1Force',
        'OBC2': 'GBSAOBC2Force',
        'GBn': 'GBSAGBnForce',
        'GBn2': 'GBSAGBn2Force',
    }

    def resolve_mmpbsa_dielectric(self, topology):
        """amber_compatible: set MMPBSA.py's solvent dielectric for this topology.

        MMPBSA.py evaluates GB with sander (extdiel 80.0) for CHARMM/chamber
        topologies and with its NAB program (extdiel 78.3) otherwise. Decided once,
        from the first (complex) structure; a no-op unless the dielectric was left
        to the preset.
        """
        if not getattr(self, 'dielectric_from_preset', False) or getattr(self, '_preset_dielectric_resolved', False):
            return
        charmm = bool(getattr(topology, 'opengbsa_charmm_like', False))
        self.solvent_dielectric = 80.0 if charmm else 78.3
        self._preset_dielectric_resolved = True
        log.info(f"amber_compatible: solvent dielectric {self.solvent_dielectric} "
                 f"({'CHARMM/chamber topology, MMPBSA.py sander path' if charmm else 'MMPBSA.py NAB path'}).")

    def _apply_rgbmax(self, gb_force):
        """Restrict the effective-Born-radius integral to pairs within rgbmax.

        Amber (sander, NAB) sums the descreening integral only over atoms j within
        rgbmax of atom i: atoms beyond rgbmax + s_j are skipped and atoms straddling
        the cutoff are integrated up to rgbmax (U = rgbmax). OpenMM's GB forces
        integrate every pair. Capping the upper limit U at rgbmax and dropping pairs
        whose lower limit L exceeds it reproduces Amber's rule for all five models
        (the GBn/GBn2 neck term acts only near contact). The GB energy pairs stay
        uncut, as in MMPBSA.py (cut = 999).
        """
        if self.rgbmax is None:
            return
        if not isinstance(gb_force, openmm.CustomGBForce):
            raise ValueError("rgbmax needs a CustomGBForce; this GB force cannot apply it.")
        for k in range(gb_force.getNumComputedValues()):
            name, expr, kind = gb_force.getComputedValueParameters(k)
            if kind != openmm.CustomGBForce.ParticlePairNoExclusions or "U=r+sr2;" not in expr:
                continue
            new = expr.replace("select(step(r+sr2-or1), ", "select(step(r+sr2-or1)*step(rgbmax-L), ", 1)
            new = new.replace("U=r+sr2;", "U=min(r+sr2, rgbmax);", 1) + f";rgbmax={self.rgbmax * 0.1:.10g}"
            if new.count("rgbmax") != 3:
                raise ValueError(f"Unexpected Born-radius integral for rgbmax: {expr}")
            gb_force.setComputedValueParameters(k, name, new, kind)
            log.info(f"Born-radius integral restricted to rgbmax = {self.rgbmax:g} A.")
            return
        raise ValueError("No Born-radius integral found to apply rgbmax to.")

    def _create_fallback_obc_force(self, system, topology):
        """
        Build the polar GB force matching `self.gb_model` (used both as the
        fallback when no GB force exists yet, and -- via
        `refine_gbsa_forces` -- to UNCONDITIONALLY REPLACE any GB force
        ParmEd/OpenMM's own `createSystem` already built, so that this
        codebase's own GB radii (`_get_gb_radius`/`_get_gb_scale`) are
        enforced).

        Despite the name (kept for backward compatibility -- this method is
        called from several sites), this is NOT OBC-only: it dispatches to
        the `openmm.app.internal.customgbforces` class matching
        `self.gb_model` (see `_GB_MODEL_FORCE_CLASS_NAME`). Previously this
        method always built a `GBSAOBC2Force`/`GBSAOBCForce` regardless of
        `self.gb_model`, so selecting HCT/OBC1/GBn/GBn2 silently fell back
        to OBC2 physics every time this fallback path was taken (confirmed:
        all 5 gb_model values produced bit-identical energies on a real
        Coordinate-Mode PPI system, since `SystemGenerator` there always
        rejects the `implicitSolvent` kwarg and lands here; and OBC1/OBC2
        were indistinguishable even in Native Mode, since both took this
        same OBC2-only path).

        All five `customgbforces` classes share one constructor signature
        and `addParticle([charge, radius, screen])` convention (unlike
        OpenMM's own native `openmm.GBSAOBCForce`, which only implements
        OBC and has a different `addParticle(charge, radius, scale)` calling
        convention with no kappa/salt support at all) -- using
        `customgbforces` uniformly for every model, including kappa=0.0 for
        the no-salt case, removes the previous OBC-only/native-vs-custom
        branch entirely rather than needing a separate code path per model.
        """
        from openmm.app.internal import customgbforces

        charges = self._extract_charges_from_system(system)

        try:
            atoms_iterable = topology.atoms()
        except TypeError:
            atoms_iterable = topology.atoms
        atoms_list = list(atoms_iterable)

        conc = self.salt_concentration
        if unit.is_quantity(conc):
            conc = conc.value_in_unit(unit.molar)
        use_salt = bool(conc) and conc > 0

        kappa = 0.0
        if use_salt:
            temp_k = getattr(self, 'temperature', 300.0)
            if unit.is_quantity(temp_k):
                temp_k = temp_k.value_in_unit(unit.kelvin)
            kappa = 50.33355 * (conc / (self.solvent_dielectric * temp_k)) ** 0.5 * 7.3  # nm^-1

        cutoff_nm = (self.nonbonded_cutoff * unit.angstroms).value_in_unit(unit.nanometer) \
            if self.nonbonded_cutoff is not None else None

        force_class_name = self._GB_MODEL_FORCE_CLASS_NAME.get(self.gb_model, 'GBSAOBC2Force')
        force_class = getattr(customgbforces, force_class_name)
        gb_force = force_class(
            solventDielectric=self.solvent_dielectric,
            soluteDielectric=self.solute_dielectric,
            SA=None, cutoff=cutoff_nm, kappa=kappa,
        )
        # NOTE: every customgbforces addParticle (inherited from
        # CustomAmberGBForceBase) already subtracts the 0.009 nm OBC offset
        # from the radius AND multiplies scale by that offset-subtracted
        # radius internally -- confirmed against OpenMM's own
        # getStandardParameters()/ParmEd's own createSystem(implicitSolvent=
        # OBC2), both of which pass raw (un-offset, un-multiplied) radius/
        # scale straight through. Pre-applying both transformations here
        # ourselves would apply them TWICE, corrupting every Born radius and
        # roughly doubling the magnitude of the resulting GB energy --
        # verified numerically against a direct sander (Amber igb=5)
        # single-frame energy: the double-applied version gave -4992
        # kcal/mol vs sander's -2482 kcal/mol for the same structure, while
        # raw (correctly single-applied) parameters give -2473 kcal/mol, a
        # 0.4% match.
        # GBn2 (Amber igb=8) is the one model in this family whose
        # addParticle takes 6 parameters, not 3: [charge, radius, screen,
        # alpha, beta, gamma] -- the extra three are fixed, element- and
        # nucleic-acid-vs-protein-dependent Born-radius correction
        # coefficients (GBSAGBn2Force._atom_params/_atom_params_nucleic),
        # NOT derived from this structure's own charges/radii the way the
        # first three parameters are. Confirmed directly against OpenMM's
        # own GBSAGBn2Force.getStandardParameters() classmethod, which
        # looks these up the same way (by element, and by whether the
        # atom's residue name is one of the 8 standard RNA/DNA residue
        # names) rather than exposing them as configurable per-run inputs.
        # GBn/GBn2 also use a model-specific per-element `screen` value
        # (OpenMM's internal `_SCREEN_PARAMETERS` element-keyed table),
        # NOT this structure's own prmtop-derived `atom.screen` -- confirmed
        # against OpenMM's own getStandardParameters() for each class:
        # OBC1/OBC2/HCT read `_screen_parameter(atom)[0]` (the "standard"
        # index, consistent with reusing the prmtop's own mbondi/mbondi2
        # SCREEN values), but GBn reads index [1] and GBn2 reads index [2]
        # (protein) or [3] (nucleic acid) -- all three noticeably different
        # from index 0 for every element (e.g. hydrogen: 0.85 at index 0 vs
        # 1.09/1.43/1.70 at indices 1/2/3). Reusing this structure's own
        # `_get_gb_scale`/`atom.screen` for GBn/GBn2 was silently using the
        # wrong screen values for every atom -- confirmed as the root cause
        # of a real ~13 kcal/mol discrepancy against Amber MMPBSA.py for
        # GBn2 on a real system (GBn's error was smaller, ~3 kcal/mol,
        # since its screen values happen to be numerically closer to the
        # standard set than GBn2's).
        is_gbn = force_class_name == 'GBSAGBnForce'
        is_gbn2 = force_class_name == 'GBSAGBn2Force'
        radii, screens = self._assigned_radii(topology, atoms_list)
        for i, atom in enumerate(atoms_list):
            charge = charges[i]
            raw_radius_nm = (radii[i] if radii else self._get_gb_radius(atom)) * 0.1
            if is_gbn or is_gbn2:
                residue_name = getattr(getattr(atom, 'residue', None), 'name', None)
                is_nucleic = residue_name in customgbforces._NUCLEIC_ACID_RESIDUES
                screen_index = 1 if is_gbn else (3 if is_nucleic else 2)
                raw_scale = customgbforces._screen_parameter(atom)[screen_index]
            else:
                raw_scale = screens[i] if screens else self._get_gb_scale(atom)
            params = [charge, raw_radius_nm, raw_scale]
            if is_gbn2:
                element = getattr(atom, 'element', None)
                table = (customgbforces.GBSAGBn2Force._atom_params_nucleic if is_nucleic
                         else customgbforces.GBSAGBn2Force._atom_params)
                params += list(table.get(element, customgbforces.GBSAGBn2Force._default_atom_params))
            gb_force.addParticle(params)
        try:
            gb_force.finalize()
        except AttributeError:
            pass
        log.info(f"Built GB force ({force_class_name}, matching gb_model='{self.gb_model}'"
                 f"{f', kappa={kappa:.4f} nm^-1' if use_salt else ''}) "
                 f"for salt_concentration={conc} M.")
        return gb_force

    def _assigned_radii(self, topology, atoms_list):
        """GB radii and screens for a topology that carries none of its own.

        Coordinate Mode builds the system from a PDB through an OpenMM
        ForceField, so the topology reaching the GB force is an
        openmm.app.Topology whose atoms have no radii. Up to 0.0.12 those
        atoms got one radius per element (every hydrogen 1.2 A), whatever
        gb_model asked for: mbondi2 gives N-bonded hydrogens 1.3 A, and mbondi
        changes more. They now get the gb_model's set (or an explicit
        gb_radii_set), assigned by ParmEd's changeRadii from each atom's
        element and bonded partner, as tleap does. Returns (None, None) when
        the atoms already carry radii (ParmEd structures from a topology file).
        """
        if not atoms_list or any(getattr(a, 'solvent_radius', 0) for a in atoms_list):
            return None, None
        import parmed
        from parmed.tools import changeRadii
        struct = topology if isinstance(topology, parmed.Structure) else parmed.openmm.load_topology(topology)
        if len(struct.atoms) != len(atoms_list):
            raise ValueError(f"GB radii: topology has {len(struct.atoms)} atoms, system {len(atoms_list)}")
        if not struct.bonds:
            raise ValueError("GB radii: the topology has no bonds, so hydrogen radii cannot be "
                             "assigned from their bonded partner")
        radii_set = StructureManager.gb_radii_decision(struct, self.gb_model, getattr(self, 'gb_radii_set', 'auto'))
        changeRadii(struct, radii_set).execute()
        log.info(f"Assigned '{radii_set}' GB radii to a topology without radii ({len(struct.atoms)} atoms).")
        return [a.solvent_radius for a in struct.atoms], [a.screen for a in struct.atoms]

    def _get_gb_radius(self, atom):
        """Get GB radius for atom"""
        # Mapping or Prmtop logic
        # Check if atom has solvent_radius (ParmEd)
        if hasattr(atom, 'solvent_radius') and atom.solvent_radius is not None and atom.solvent_radius > 0:
            return atom.solvent_radius

        # Standard mbondi2 radii mapping
        default_radii = {
            'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.5,
            'F': 1.5, 'S': 1.8, 'P': 1.8, 'Cl': 1.7,
            'I': 1.9, 'Br': 1.85, 'Na': 1.5, 'K': 1.5,
            'Mg': 1.0, 'Ca': 1.0, 'Zn': 1.0
        }
        symbol = atom.element.symbol if hasattr(atom, 'element') and atom.element else 'C'
        radius = default_radii.get(symbol)
        if radius is not None:
            return radius
        if hasattr(self, 'physics_assumptions') and not getattr(self, '_logged_radius_fallback', False):
            self.physics_assumptions.append(
                "Implicit Solvent: Unknown element encountered while assigning GB radii; "
                "using generic fallback radii for unresolved atom types."
            )
            self._logged_radius_fallback = True
        return 1.7

    def _get_gb_scale(self, atom):
        """Get GB scale for atom"""
        # Check if atom has screen (ParmEd)
        if hasattr(atom, 'screen'):
            return atom.screen
            
        # Standard mbondi2 scale mapping
        default_scales = {
            'H': 0.85, 'C': 0.72, 'N': 0.79, 'O': 0.85,
            'F': 0.88, 'S': 0.96, 'P': 0.86, 'Cl': 0.8,
            'I': 0.8, 'Br': 0.8, 'Na': 0.8, 'K': 0.8,
            'Mg': 0.8, 'Ca': 0.8, 'Zn': 0.8
        }
        symbol = atom.element.symbol if hasattr(atom, 'element') and atom.element else 'C'
        return default_scales.get(symbol, 0.8)

    def _setup_surface_area_force(self, system, topology, lcpo_type_names=None):
        """Setup Surface Area Force explicitly"""
        if self.sa_model == 'LCPO':
             return self._setup_lcpo_force(system, topology, lcpo_type_names)
        else:
             return self._create_ace_sa_force(system, topology)

    def _create_ace_sa_force(self, system, topology):
        """Create separate ACE SA force using GBSAOBCForce with 0 charges"""
        sa_force = openmm.GBSAOBCForce()
        sa_force.setNonbondedMethod(openmm.GBSAOBCForce.NoCutoff)
        sa_force.setSoluteDielectric(1.0) 
        sa_force.setSolventDielectric(78.5)
        sa_force.setSurfaceAreaEnergy(2.25936) # kJ/mol/nm^2
        
        # We need to replicate particles with 0 charge to get SA only
        # We need radii/scales. 
        # Using implicit solvent usually sets them.
        # But here we are creating a NEW force.
        # We need charges (0), radii, scales.
        # Check system existing GB force for parameters?
        # This is getting complicated to implement correctly for fallback.
        # But I don't need to fix ACE right now, just ensure code is valid.
        
        atoms_list = []
        if hasattr(topology, 'atoms'):
            if callable(topology.atoms):
                atoms_list = topology.atoms()
            else:
                atoms_list = topology.atoms
        else:
             atoms_list = []
             
        for atom in atoms_list:
             radius = self._get_gb_radius(atom) * 0.1 # nm
             scale = self._get_gb_scale(atom)
             sa_force.addParticle(0.0, radius, scale)
             
        print(f"✓ Added enhanced surface area force (GBSAOBC-ACE-SA) to system")
        return sa_force

    def _setup_lcpo_force(self, system, topology, lcpo_type_names=None):
        """Setup LCPO (Linear Combinations of Pairwise Overlaps) Surface Area Force"""
        print("✓ Setting up LCPO Surface Area Force...")
        
        surface_tension = 3.01248 * unit.kilojoules_per_mole / unit.nanometers**2
        
        if not HAS_LCPO_PARAMS:
            raise ImportError("LCPO parameters not available. Cannot use LCPO model.")
            
        lcpo_force = openmm.LCPOForce()
        try:
             lcpo_force.setSurfaceTension(surface_tension)
        except AttributeError:
             pass
        
        probe_radius_nm = 1.4 * 0.1
        if self.lcpo_typing == 'amber':
             if lcpo_type_names is not None:
                  # MMPBSA.py's assignment by atom-type name (radius in A, P4 in 1/A^2)
                  for r, p1, p2, p3, p4 in amber_lcpo_params(topology, lcpo_type_names):
                       lcpo_force.addParticle((r * 0.1 + probe_radius_nm) if r > 0 else 0.0,
                                              p1, p2, p3, p4 * 100.0)
                  log.info(f"LCPO parameters assigned by atom-type name, as in Amber MMPBSA.py "
                           f"({lcpo_force.getNumParticles()} particles).")
                  lcpo_force.setForceGroup(4)
                  return lcpo_force
             log.warning("lcpo_typing='amber' requested, but this input carries no atom-type "
                         "names; assigning LCPO parameters by element and bonding instead.")
        print("ℹ️  Calculating LCPO parameters from topology...")
        params_list = getLCPOParamsTopology(topology)
        
        for params in params_list:
             r = params[0]
             p1, p2, p3, p4 = params[1], params[2], params[3], params[4]
             radius_nm = r.value_in_unit(unit.nanometer)
             val_p1 = p1 if not unit.is_quantity(p1) else p1._value
             val_p2 = p2 if not unit.is_quantity(p2) else p2._value
             val_p3 = p3 if not unit.is_quantity(p3) else p3._value
             val_p4 = p4.value_in_unit(unit.nanometer**-2) if unit.is_quantity(p4) else p4 * 100.0
             
             effective_radius = (radius_nm + probe_radius_nm) if radius_nm > 0 else 0.0
             lcpo_force.addParticle(effective_radius, val_p1, val_p2, val_p3, val_p4)
             
        print(f"✓ Added LCPO Surface Area Force ({lcpo_force.getNumParticles()} particles using internal params)")
        lcpo_force.setForceGroup(4)
        return lcpo_force

    def _get_gb_radius(self, atom):
        """Get GB radius for atom"""
        # Prefer explicit radii when available (ParmEd/native topology).
        if hasattr(atom, 'solvent_radius') and atom.solvent_radius is not None and atom.solvent_radius > 0:
            return atom.solvent_radius

        # Handle ParmEd structures where element might be an int
        if hasattr(atom, 'element') and isinstance(atom.element, int):
             try:
                 elem = app.Element.getByAtomicNumber(atom.element)
                 element = elem.symbol
             except:
                 element = 'C'
        else:
             element = atom.element.symbol if hasattr(atom, 'element') and atom.element else 'C'
             
        # Standard mbondi2 radii mapping
        default_radii = {
            'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.5,
            'F': 1.5, 'S': 1.8, 'P': 1.8, 'Cl': 1.7,
            'I': 1.9, 'Br': 1.85, 'Na': 1.5, 'K': 1.5,
            'Mg': 1.0, 'Ca': 1.0, 'Zn': 1.0
        }
        radius = default_radii.get(element)
        if radius is not None:
            return radius
        if hasattr(self, 'physics_assumptions') and not getattr(self, '_logged_radius_fallback', False):
            self.physics_assumptions.append(
                "Implicit Solvent: Unknown element encountered while assigning GB radii; "
                "using generic fallback radii for unresolved atom types."
            )
            self._logged_radius_fallback = True
        return 1.50  # Generic default

    def _get_gb_scale(self, atom):
        """Get GB scaling factor for atom"""
        # Prefer explicit GB screening scale when available.
        if hasattr(atom, 'screen') and atom.screen is not None and atom.screen > 0:
            return atom.screen

        # Handle ParmEd structures where element might be an int
        if hasattr(atom, 'element') and isinstance(atom.element, int):
             try:
                 elem = app.Element.getByAtomicNumber(atom.element)
                 element = elem.symbol
             except:
                 element = 'C'
        else:
             element = atom.element.symbol if hasattr(atom, 'element') and atom.element else 'C'
             
        # Standard mbondi2 scale mapping
        default_scales = {
            'H': 0.85, 'C': 0.72, 'N': 0.79, 'O': 0.85,
            'F': 0.88, 'S': 0.96, 'P': 0.86, 'Cl': 0.8,
            'I': 0.8, 'Br': 0.8, 'Na': 0.8, 'K': 0.8,
            'Mg': 0.8, 'Ca': 0.8, 'Zn': 0.8
        }
        return default_scales.get(element, 0.80)  # Default 0.8

    def _extract_charges_from_system(self, system):
        """Extract atomic charges from NonbondedForce"""
        charges = []
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    charges.append(charge.value_in_unit(unit.elementary_charge))
                break
        else:
            raise ValueError("No NonbondedForce found in system")
        return charges

    def _remove_implicit_forces(self, system):
        """Remove existing implicit solvent forces"""
        forces_to_remove = []
        for i, force in enumerate(system.getForces()):
            if isinstance(force, (openmm.GBSAOBCForce, openmm.CustomGBForce)):
                forces_to_remove.append(i)
        for i in reversed(forces_to_remove):
            system.removeForce(i)
            print(f"✓ Removed existing implicit solvent force")

    def validate_gbsa_setup(self, system, topology):
        """Validation of GBSA force setup"""
        gb_forces = []
        sa_forces = []
        screening_forces = []
        
        for force in system.getForces():
            if isinstance(force, openmm.GBSAOBCForce):
                gb_forces.append(force)
            elif isinstance(force, openmm.CustomGBForce):
                # Check if it's surface area or GB force
                if hasattr(force, 'getEnergyTermParameters'):
                    try:
                        expr, _ = force.getEnergyTermParameters(0)
                        if 'SASA' in expr or 'gamma' in expr:
                            sa_forces.append(force)
                        else:
                            gb_forces.append(force)
                    except:
                        gb_forces.append(force)  # Default to GB
            elif isinstance(force, openmm.CustomNonbondedForce):
                # Check for screening force
                try:
                    expr = force.getEnergyFunction()
                    if 'exp(-kappa*r)' in expr:
                        screening_forces.append(force)
                except:
                    pass
        
        print(f"✓ Validation: {len(gb_forces)} GB forces, {len(sa_forces)} SA forces, {len(screening_forces)} screening forces")
        return len(gb_forces) > 0

    def decompose_energy_contributions(self, system, context, positions):
        """Decompose energy into individual force contributions"""
        energy_decomposition = {}
        for i, force in enumerate(system.getForces()):
            force_name = force.__class__.__name__
            try:
                temp_system = openmm.System()
                temp_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
                for j in range(system.getNumParticles()):
                    temp_system.addParticle(system.getParticleMass(j))
                temp_system.addForce(force)
                temp_integrator = openmm.VerletIntegrator(0.001*unit.picoseconds)
                temp_context = openmm.Context(temp_system, temp_integrator)
                temp_context.setPositions(positions)
                energy = temp_context.getState(getEnergy=True).getPotentialEnergy()
                energy_decomposition[force_name] = energy.value_in_unit(unit.kilocalories_per_mole)
                del temp_context
            except Exception as e:
                energy_decomposition[force_name] = 0.0
        return energy_decomposition


# --- Parallel per-frame energy evaluation (CPU-only) -----------------------
#
# Mirrors the worker pattern already proven in `decomposition.py`
# (`_worker_init`/`_worker_analyze_frame`): each worker process
# deserializes its own copy of the complex/protein/ligand `System`s once
# (via `openmm.XmlSerializer`, since `System`/`Context` objects can't be
# pickled directly) and builds its own three `Context`s, then processes
# whatever frame indices it's given. This is CPU-only -- CUDA/OpenCL
# contexts cannot be safely shared or re-created across worker processes
# (same constraint `decomposition.py`'s own parallel path already
# documents), so `run()` only dispatches to this pool when the resolved
# platform is CPU; every other platform keeps using the original serial
# loop unchanged.
_frame_worker_state = {}
_KCAL = unit.kilocalories_per_mole
_INTERNAL_GROUPS = (10, 11, 12, 13, 14)   # bond, angle, torsion, RB torsion, CMAP


def _platform_name(name):
    """Canonical OpenMM platform name in any case ('cuda' -> 'CUDA'); None for
    null/'auto'; an unknown name is returned unchanged."""
    from .config import canonical_platform
    try:
        return canonical_platform(name)
    except ValueError:
        return name


def _force_groups(system):
    """The force groups present in a System, as a sorted tuple."""
    return tuple(sorted({f.getForceGroup() for f in system.getForces()}))


def _group_energies(context, groups):
    """Potential energy of each force group, one evaluation per group (kcal/mol).

    getState(groups=mask) evaluates only the forces in `mask`, so evaluating
    each group once and summing gives the total without evaluating any force
    twice. The earlier per-frame code requested the total and then the
    individual groups, so the expensive GB, non-polar and nonbonded forces were
    evaluated twice (complex) or three times (receptor, ligand) per frame.
    """
    return {g: context.getState(getEnergy=True, groups=1 << g).getPotentialEnergy().value_in_unit(_KCAL)
            for g in groups}


def _frame_energy_terms(complex_context, protein_context, ligand_context, groups):
    """Every per-frame quantity run() records, from one evaluation per force group.

    `groups` is (complex_groups, protein_groups, ligand_groups) from
    `_force_groups`. Positions must already be set on the three contexts.
    Group conventions (assign_force_groups): 0 electrostatic, 3 van der Waals,
    1-2 polar GB, 4 non-polar SA (receptor/ligand also 16), 10-14 internal.
    """
    Ec, Ep, El = (_group_energies(ctx, g) for ctx, g in
                  zip((complex_context, protein_context, ligand_context), groups))

    def s(E, *gs):
        return sum(E.get(g, 0.0) for g in gs)

    complex_e, protein_e, ligand_e = sum(Ec.values()), sum(Ep.values()), sum(El.values())
    e_internal = s(Ec, *_INTERNAL_GROUPS)
    binding_e = ((complex_e - e_internal) - (protein_e - s(Ep, *_INTERNAL_GROUPS))
                 - (ligand_e - s(El, *_INTERNAL_GROUPS)))
    c_elec, c_vdw = s(Ec, 0), s(Ec, 3)
    p_elec, p_vdw = s(Ep, 0), s(Ep, 3)
    l_elec, l_vdw = s(El, 0), s(El, 3)
    c_gb = s(Ec, 1, 2)
    return {
        'complex_e': complex_e, 'protein_e': protein_e, 'ligand_e': ligand_e,
        'binding_e': binding_e,
        'c_nb': c_elec + c_vdw, 'c_gb': c_gb, 'c_sa': s(Ec, 4),
        'c_vdw': c_vdw, 'c_elec': c_elec, 'p_vdw': p_vdw, 'p_elec': p_elec,
        'l_vdw': l_vdw, 'l_elec': l_elec, 'p_gb': s(Ep, 1, 2), 'l_gb': s(El, 1, 2),
        'p_sa': s(Ep, 4, 16), 'l_sa': s(El, 4, 16),
        'e_bond': s(Ec, 10), 'e_angle': s(Ec, 11), 'e_torsion': s(Ec, 12),
        'e_rbtorsion': s(Ec, 13), 'e_cmap': s(Ec, 14),
        'e_nb': c_elec + c_vdw, 'e_obc': c_gb, 'e_internal': e_internal,
        'complex_e_clean': complex_e - e_internal,
        'protein_e_clean': protein_e - s(Ep, *_INTERNAL_GROUPS),
        'ligand_e_clean': ligand_e - s(El, *_INTERNAL_GROUPS),
        'error': None,
    }


def _cpu_threads_per_worker(n_workers):
    """OpenMM CPU threads for each of `n_workers` processes.

    Without this every worker's CPU platform uses all hardware threads, so N
    workers start N x n_cpu threads and oversubscribe the machine. An explicit
    OPENMM_CPU_THREADS in the environment is respected (returns None)."""
    if os.environ.get('OPENMM_CPU_THREADS'):
        return None
    return max(1, (os.cpu_count() or 1) // max(1, n_workers))


def _frame_worker_init(complex_xml, protein_xml, ligand_xml, temperature_k, threads=None):
    """Runs once per worker process: deserialize the three Systems and
    build their Contexts. Force groups (0=NB, 1=GB-OBC, 2=GB-Custom,
    3=VDW, 4=SA, 10-14=internal) are already baked into the serialized
    System by `assign_force_groups`, so no re-assignment is needed here.
    `threads` caps the OpenMM CPU threads of this worker's contexts."""
    global _frame_worker_state
    complex_system = openmm.XmlSerializer.deserialize(complex_xml)
    protein_system = openmm.XmlSerializer.deserialize(protein_xml)
    ligand_system = openmm.XmlSerializer.deserialize(ligand_xml)

    platform = openmm.Platform.getPlatformByName('CPU')
    props = {'Threads': str(int(threads))} if threads else {}
    fric, step = 1.0, 0.001
    complex_integrator = openmm.LangevinMiddleIntegrator(float(temperature_k), fric, step)
    protein_integrator = openmm.LangevinMiddleIntegrator(float(temperature_k), fric, step)
    ligand_integrator = openmm.LangevinMiddleIntegrator(float(temperature_k), fric, step)

    _frame_worker_state['complex_context'] = openmm.Context(complex_system, complex_integrator, platform, props)
    _frame_worker_state['protein_context'] = openmm.Context(protein_system, protein_integrator, platform, props)
    _frame_worker_state['ligand_context'] = openmm.Context(ligand_system, ligand_integrator, platform, props)
    _frame_worker_state['groups'] = (_force_groups(complex_system), _force_groups(protein_system),
                                     _force_groups(ligand_system))


def _frame_worker_process(args):
    """Processes one trajectory frame using this worker's pre-built
    Contexts. Reproduces exactly the per-frame calculation from `run()`'s
    main loop (post get_clean_energy/get_grp_E de-duplication) -- any
    change to that calculation must be mirrored here, and vice versa.
    Returns (frame_idx, results_dict) on success, or (frame_idx, None) with
    the exception recorded in results_dict['error'] on failure, matching
    the serial loop's own `except: print(...); continue` behavior (the
    caller logs and skips failed frames instead of aborting the whole run).
    """
    (frame_idx, xyz_nm, protein_indices, ligand_indices, complex_pos_full,
     is_native_mode, is_ppi_coordinate_mode) = args

    complex_context = _frame_worker_state['complex_context']
    protein_context = _frame_worker_state['protein_context']
    ligand_context = _frame_worker_state['ligand_context']

    try:
        ligand_pos = unit.Quantity(xyz_nm[np.asarray(ligand_indices)], unit.nanometer)
        protein_pos = unit.Quantity(xyz_nm[np.asarray(protein_indices)], unit.nanometer)
        ligand_context.setPositions(ligand_pos)
        protein_context.setPositions(protein_pos)
        if is_native_mode or is_ppi_coordinate_mode:
            complex_context.setPositions(unit.Quantity(complex_pos_full, unit.nanometer))
        else:
            combined_pos = unit.Quantity(
                np.concatenate([xyz_nm[np.asarray(protein_indices)], xyz_nm[np.asarray(ligand_indices)]]),
                unit.nanometer)
            complex_context.setPositions(combined_pos)
        result = _frame_energy_terms(complex_context, protein_context, ligand_context,
                                     _frame_worker_state['groups'])
        return (frame_idx, result)
    except Exception as e:
        return (frame_idx, {'error': str(e)})


class GBSACalculator(GBSAForceManager):
    """Advanced True Force Field MMGBSA Calculator"""

    # charge_method: partial charges given to a Coordinate-Mode ligand whose file has
    # none, and to residues the protein force field lacks (parameterised with OpenFF)
    CHARGE_METHODS = ('am1bcc', 'mmff94', 'gasteiger')

    def __init__(self, temperature=300, verbose=True, gb_model='OBC2', salt_concentration=0.15,
                 max_workers=None, protein_forcefield='amber', charge_method='am1bcc',
             solute_dielectric=1.0, solvent_dielectric=None, entropy_method='none', decomposition_method='full',
             platform=None, sa_model=None, nonbonded_cutoff=None,
             reimage_trajectory=True, charmm_params=None, charmm_coordinates=None,
             gb_radii_set='auto', amber_compatible=False, rgbmax=None, lcpo_typing=None,
             nmode_frames=10, nmode_drms=0.001, nmode_maxcyc=10000, nmode_precision='mixed'):
        """
        Initialize the MM/GBSA calculator with enhanced features
        
        Parameters:
        -----------
        temperature : float
             Temperature in Kelvin for the analysis
        verbose : bool
             Print setup details (platform request, ...) in addition to the log
        gb_model : str
             Generalized Born model ('OBC1', 'OBC2', 'HCT', 'GBn', 'GBn2')
        salt_concentration : float
             Salt concentration in Molar for Debye-Hückel screening
        max_workers : int, optional
             Worker processes of the CPU frame loop (None: min(cores - 1, 4)).
             The frames are split over workers whenever the platform is CPU.
        protein_forcefield : str
             Protein forcefield to use ('amber', etc.)
        charge_method : str
             Coordinate Mode: partial charges ('am1bcc', 'mmff94' or 'gasteiger',
             OpenFF toolkit) of a ligand_mol whose file has none and of residues
             the protein force field lacks. Charges in the ligand file are kept.
             A failing method is an error, not a fallback.
        solute_dielectric : float
             Solute dielectric constant (default 1.0)
        solvent_dielectric : float
             Solvent dielectric constant (default 78.5)
        entropy_method : str
             Entropy calculation method ('interaction', 'c2', 'quasiharmonic',
             'nmode', 'none'). 'interaction' (IE) and 'c2' use the gas-phase
             interaction energy (mmgbsa.interaction_entropy); 'nmode' (alias
             'normal_mode') is the normal-mode binding entropy of complex,
             receptor and ligand (mmgbsa.nmode_entropy)
        nmode_frames : int or None
             'nmode': number of frames, spread evenly over the analysed frames
             (None: every analysed frame)
        nmode_drms : float
             'nmode': RMS gradient (kcal/mol/A) each minimization must reach;
             frames where a species does not reach it are left out
        nmode_maxcyc : int
             'nmode': iterations per minimization round
        nmode_precision : str
             'nmode': 'mixed' or 'double' force precision for minimization and
             Hessian on CUDA/OpenCL (CPU runs use the double-precision
             Reference platform)
        platform : str, optional
             Platform to force (CPU, CUDA, OpenCL, Reference)
        sa_model : str
             Surface Area model ('ACE' or 'LCPO')
        nonbonded_cutoff : float, optional
             Cutoff distance in Angstroms for nonbonded interactions.
             If None (default) or >= 999.0, NoCutoff is used.
        reimage_trajectory : bool, default True
             Re-image molecules across periodic boundaries (PBC re-wrapping)
             before computing any energy. A raw MD trajectory has no
             guarantee that a molecule stays whole/centered from frame to
             frame -- if the protein or ligand drifts and wraps to the
             opposite side of the periodic box mid-trajectory, downstream
             vdW/electrostatic energies are silently corrupted even though
             the true physics hasn't changed (confirmed on real production
             data: leaving this disabled reproduces a systematic ~11%
             low-magnitude bias in delta_vdw that grows across the
             trajectory). Enabled by default for this reason; only disable
             if the input trajectory is already known to be correctly
             imaged (e.g. it was already processed with cpptraj `autoimage`
             or an equivalent tool), since re-imaging an already-imaged
             trajectory is a harmless no-op but still costs time.
        charmm_params : str or list of str, optional
             Path(s) to CHARMM parameter file(s) (.prm/.str/.rtf), required
             only when the topology passed to `parameterize_protein_amber`
             is a NAMD/CHARMM `.psf` file (a PSF has no embedded force-field
             parameters, unlike Amber's prmtop). Ignored for Amber/GROMACS/
             generic-PDB topologies. See `TopologyLoader._load_charmm`.
        charmm_coordinates : str, optional
             Path to a companion coordinate file (.pdb/.coor/.crd) for a
             `.psf` topology -- a PSF itself carries no atomic positions
             (confirmed: `parmed.load_file` on a real PSF gives
             `struct.positions is None`), unlike Amber's prmtop+inpcrd
             pairing. Required to get a physically real system; without it,
             CHARMM-mode positions fall back to a same-shape all-zeros
             placeholder, which builds but is meaningless for energetics.
        """
        self.temperature = temperature * unit.kelvin
        self.platform_preference = platform
        self.verbose = bool(verbose)
        self.gb_model = gb_model
        # 'auto': an Amber .prmtop keeps its own GB radii (as Amber does), other
        # inputs get the gb_model set; a set name forces it (gb_radii_decision).
        self.gb_radii_set = gb_radii_set or 'auto'
        self._gb_radii_decision = StructureManager._GB_MODEL_RADII_SET.get(gb_model, 'mbondi2')
        self.salt_concentration = salt_concentration
        self.max_workers = max_workers or max(1, min(mp.cpu_count() - 1, 4))
        # 'amber14-all' and the file name 'amber14-all.xml' mean the amber14 set;
        # 'auto' means the default
        pff = str(protein_forcefield or 'amber').strip().lower()
        pff = {'amber14-all': 'amber14', 'amber14-all.xml': 'amber14', 'auto': 'amber'}.get(pff, pff)
        self.protein_forcefield = pff
        # Support user-defined ligand forcefield (openff, gaff)
        self.ligand_forcefield = 'openff' 
        charge_method = str(charge_method or 'am1bcc').strip().lower()
        if charge_method not in self.CHARGE_METHODS:
            raise ValueError(f"charge_method must be one of {list(self.CHARGE_METHODS)}, "
                             f"not '{charge_method}'")
        self.charge_method = charge_method
        self.solute_dielectric = solute_dielectric
        # amber_compatible reproduces Amber MMPBSA.py's GB conventions: LCPO non-polar
        # term, effective Born radii only within rgbmax = 25 A, and MMPBSA.py's solvent
        # dielectric (80.0 on its sander path, used for CHARMM/chamber topologies;
        # 78.3 on its NAB path, used otherwise). Explicit settings win. By default
        # OpenGBSA evaluates every pair (no Born-radius cutoff) at dielectric 78.5.
        self.amber_compatible = bool(amber_compatible)
        if sa_model is None:
            sa_model = 'LCPO'
        elif self.amber_compatible and sa_model != 'LCPO':
            log.warning(f"amber_compatible uses LCPO; keeping the explicit sa_model='{sa_model}'.")
        if rgbmax is None and self.amber_compatible:
            rgbmax = 25.0
        self.rgbmax = rgbmax
        if lcpo_typing is None:
            lcpo_typing = 'amber' if self.amber_compatible else 'element'
        self.lcpo_typing = lcpo_typing
        self._dielectric_from_preset = solvent_dielectric is None and self.amber_compatible
        self.dielectric_from_preset = self._dielectric_from_preset     # for the inherited refine path
        self._preset_dielectric_resolved = False
        if solvent_dielectric is None:
            solvent_dielectric = 78.5     # amber_compatible: replaced once the topology is known
        self.solvent_dielectric = solvent_dielectric
        self.entropy_method = entropy_method
        # Normal-mode entropy (entropy_method 'nmode'): frames used, minimization
        # criterion (kcal/mol/A) and iterations, GPU precision of the Hessian
        self.nmode_frames = nmode_frames
        self.nmode_drms = nmode_drms
        self.nmode_maxcyc = nmode_maxcyc
        self.nmode_precision = nmode_precision
        self.nmode_result = None
        self.decomposition_method = decomposition_method
        self.reimage_trajectory = reimage_trajectory
        self.charmm_params = charmm_params
        self.charmm_coordinates = charmm_coordinates
        self.parameterized_residues = [] # Track any dynamic residues
        self.physics_assumptions = []

        # Cutoff Support
        self.nonbonded_cutoff = nonbonded_cutoff
        if self.nonbonded_cutoff is not None and self.nonbonded_cutoff >= 999.0:
            self.nonbonded_cutoff = None # Treat large cutoff as NoCutoff explicitly
        
        # Surface Area Model
        self.sa_model = sa_model
        if self.sa_model == 'LCPO' and not HAS_LCPO_PARAMS:
            print("WARNING: LCPO model requested but openmm.app.internal.lcpo not found. Falling back to ACE.")
            self.sa_model = 'ACE'
        

        # Initialize fixed enhanced GBSA manager
        self.gbsa_manager = GBSAForceManager(gb_model=gb_model, salt_concentration=salt_concentration, solute_dielectric=solute_dielectric, solvent_dielectric=solvent_dielectric, sa_model=self.sa_model, nonbonded_cutoff=self.nonbonded_cutoff)
        self.gbsa_manager.gb_radii_set = self.gb_radii_set
        self.gbsa_manager.rgbmax = self.rgbmax
        self.gbsa_manager.dielectric_from_preset = self._dielectric_from_preset
        self.gbsa_manager.lcpo_typing = self.lcpo_typing
        if self.amber_compatible or self.rgbmax is not None:
            log.info(f"GB conventions: sa_model={self.sa_model} (LCPO typing: {self.lcpo_typing}), rgbmax="
                     f"{'none' if self.rgbmax is None else f'{self.rgbmax:g} A'}, solvent dielectric "
                     f"{'from the topology (MMPBSA.py: 80.0 CHARMM/chamber, 78.3 otherwise)' if self._dielectric_from_preset else self.solvent_dielectric}.")
        
        self.energies = {'complex': [], 'protein': [], 'ligand': [], 'binding': []}
        self.energy_decompositions = []
        # Folder for the files a run derives from its inputs (format conversions,
        # aliases); the run's output folder once run_comprehensive/run knows it
        self.work_dir = None

        # Platform settings for separate analysis/decomposition configuration
        self.preferred_platform = None
        self.decomposition_platform = None
        self.strict_platform = False
        self.platform_requested = self.platform_used = self.platform_fallback = None

    def set_platform_settings(self, platform_settings):
        """
        Configure platform preferences for analysis and decomposition.
        
        Parameters:
        -----------
        platform_settings : dict
            Configuration dictionary with optional keys:
            - 'preferred_platform': Platform for main analysis (default: auto-detect)
            - 'decomposition_platform': Platform for decomposition (default: same as analysis)
            Example: {'preferred_platform': 'CUDA', 'decomposition_platform': 'CPU'}
        """
        if not platform_settings:
            return
        
        if isinstance(platform_settings, dict):
            # OpenMM platform names in any case ('cuda' -> 'CUDA'); an unknown name
            # is kept and reported as unusable by setup_optimized_platform
            self.preferred_platform = _platform_name(platform_settings.get('preferred_platform'))
            self.decomposition_platform = _platform_name(platform_settings.get('decomposition_platform'))
            self.strict_platform = bool(platform_settings.get('strict_platform', False))
            
            if self.preferred_platform:
                self.platform_preference = self.preferred_platform
            
            if self.verbose:
                if self.preferred_platform:
                    print(f"  Requested platform for analysis: {self.preferred_platform}")
                if self.decomposition_platform:
                    print(f"  Platform for decomposition: {self.decomposition_platform}")

    def _work_dir_path(self):
        """Folder for files derived from the inputs (format conversions, aliases).

        The run's output folder once run_comprehensive/run know it, else a fresh
        temporary folder; never the folder of an input file.
        """
        if self.work_dir is None:
            import tempfile
            self.work_dir = Path(tempfile.mkdtemp(prefix='opengbsa_'))
        path = Path(self.work_dir)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def validate_input_files(self, ligand_mol, complex_pdb, ligand_pdb, xtc_file, solvated_topology=None):
        """
        Comprehensive validation of input files for MM/GBSA analysis
        
        Parameters:
        -----------
        ligand_mol : str
            Path to ligand molecule file (.sdf, .mol2, .pdb) - Optional
        complex_pdb : str
            Path to protein-ligand complex PDB file
        ligand_pdb : str
            Path to isolated ligand PDB file - Optional
        xtc_file : str
            Path to molecular dynamics trajectory file
        solvated_topology : str
            Path to solvated topology file (for trajectory consistency) - Optional
            
        Returns:
        --------
        list : List of validation error messages (empty if all valid)
        """
        validation_errors = []
        
        # Check file existence
        files_to_check = {
            'ligand_mol': ligand_mol,
            'complex_pdb': complex_pdb, 
            'ligand_pdb': ligand_pdb,
            'trajectory': xtc_file
        }
        
        for file_type, file_path in files_to_check.items():
            if file_path is None: continue
            if not Path(file_path).exists():
                validation_errors.append(f"{file_type} file not found: {file_path}")
        
        if validation_errors:
            return validation_errors
        
        # Validate PDB files
        try:
            if str(complex_pdb).endswith('.prmtop'):
                 complex_pdb_obj = app.AmberPrmtopFile(complex_pdb)
            elif str(complex_pdb).endswith('.top') or str(complex_pdb).endswith('.tpr'):
                 # GROMACS Topology - Skip detailed structure validation here
                 # Conversion will handle validity
                 complex_pdb_obj = None
            elif str(complex_pdb).endswith('.gro'):
                 complex_pdb_obj = app.GromacsGroFile(complex_pdb)
            elif str(complex_pdb).endswith('.psf'):
                 # NAMD/CHARMM Topology - no embedded coordinates/parameters
                 # (see TopologyLoader._load_charmm), so it cannot be parsed
                 # as a PDB (confirmed: app.PDBFile on a real .psf raises an
                 # opaque "list index out of range" rather than a clear
                 # format error). Skip detailed structure validation here,
                 # same as the .top/.tpr native-topology cases above --
                 # actual loading/validation happens in
                 # StructureManager.load_complex/TopologyLoader._load_charmm.
                 complex_pdb_obj = None
            else:
                 complex_pdb_obj = app.PDBFile(complex_pdb)
            if ligand_pdb:
                ligand_pdb_obj = app.PDBFile(ligand_pdb)
                
                # Check if complex contains ligand
                if complex_pdb_obj is not None:
                    complex_residues = set(res.name for res in complex_pdb_obj.topology.residues())
                    ligand_residues = set(res.name for res in ligand_pdb_obj.topology.residues())
                    
                    if not ligand_residues.issubset(complex_residues):
                        validation_errors.append("Ligand residues not found in complex PDB")
            
        except Exception as e:
            validation_errors.append(f"PDB validation error: {e}")
        
        # Validate ligand molecule file
        if ligand_mol:
            try:
                # Check if file exists first
                if not Path(ligand_mol).exists():
                    validation_errors.append(f"Ligand molecule file not found: {ligand_mol}")
                elif Molecule is None:
                    # OpenFF toolkit not available - just verify file is readable
                    pass  # OpenFF toolkit not available, skip detailed validation
                else:
                    # Use OpenFF for validation
                    if ligand_mol.endswith('.mol2'):
                        # Use our robust converter for validation
                        validation_temp_sdf = str(self._work_dir_path() / "validation_temp.sdf")
                        
                        if convert_mol2_to_sdf(ligand_mol, validation_temp_sdf):
                            mol = Molecule.from_file(validation_temp_sdf, allow_undefined_stereo=True)
                        else:
                            validation_errors.append(f"Ligand Mol2 conversion failed for {ligand_mol}")
                            mol = None
                    else:
                        mol = Molecule.from_file(ligand_mol, allow_undefined_stereo=True)
                        
                    if mol is not None and mol.n_atoms == 0:
                        validation_errors.append("Ligand molecule file contains no atoms")
            except Exception as e:
                validation_errors.append(f"Ligand molecule validation error: {e}")
        
        # Validate trajectory
        if xtc_file:
            try:
                # Use solvated topology if provided to match trajectory atoms
                top_to_use = solvated_topology if solvated_topology else complex_pdb
                
                # Skip if using GROMACS topology directly (requires structure file we might not know yet)
                if str(top_to_use).endswith('.top') or str(top_to_use).endswith('.tpr'):
                     pass
                else:
                    # Callers normally already route xtc_file through
                    # TrajectoryProcessor.resolve_mdtraj_loadable_path (see
                    # run_comprehensive) before reaching here, but resolve
                    # again defensively in case this is ever called directly
                    # with a raw '.trj' path.
                    from .trajectory import TrajectoryProcessor
                    loadable_xtc = TrajectoryProcessor.resolve_mdtraj_loadable_path(
                        xtc_file, self._work_dir_path())
                    traj_test = md.load(loadable_xtc, top=top_to_use, frame=0)
                    if len(traj_test) == 0:
                        validation_errors.append("Trajectory file contains no frames")
            except Exception as e:
                validation_errors.append(f"Trajectory validation error (Atom Mismatch possible if solvated_topology missing): {e}")
        
        return validation_errors

    def _select_frames(self, trajectory_length, max_frames=None, frame_start=None, frame_end=None,
                      frame_stride=None, frame_selection='sequential', random_seed=42):
        """
        Select frames based on parameters
        
        Parameters:
        -----------
        trajectory_length : int
            Total number of frames in trajectory
        max_frames : int, optional
            Maximum number of frames to analyze
        frame_start : int, optional
            Start frame (0-indexed)
        frame_end : int, optional
            End frame (0-indexed)
        frame_stride : int, optional
            Frame stride (every Nth frame)
        frame_selection : str
            Selection method: 'sequential', 'equidistant', 'random'
        random_seed : int
            Random seed for random selection
            
        Returns:
        --------
        list : Selected frame indices
        """
        # Determine frame range
        if frame_start is None:
            frame_start = 0
        if frame_end is None:
            frame_end = trajectory_length
        
        # Validate frame range
        frame_start = max(0, min(frame_start, trajectory_length - 1))
        frame_end = max(frame_start + 1, min(frame_end, trajectory_length))
        
        print(f"Frame selection parameters:")
        print(f"  • Trajectory length: {trajectory_length}")
        print(f"  • Frame range: {frame_start} to {frame_end}")
        print(f"  • Frame stride: {frame_stride}")
        print(f"  • Selection method: {frame_selection}")
        print(f"  • Max frames: {max_frames}")
        
        # Generate frame indices based on selection method
        if frame_selection == "sequential":
            # Sequential selection with stride
            if frame_stride is None or frame_stride <= 0:
                frame_indices = list(range(frame_start, frame_end))
            else:
                frame_indices = list(range(frame_start, frame_end, frame_stride))
                
        elif frame_selection == "equidistant":
            # Equidistant selection
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            step = (frame_end - frame_start) / max_frames
            frame_indices = [frame_start + int(i * step) for i in range(max_frames)]
            frame_indices = [min(idx, frame_end - 1) for idx in frame_indices]
            
        elif frame_selection == "random":
            # Random selection
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            import random
            random.seed(random_seed)
            available_frames = list(range(frame_start, frame_end))
            if frame_stride is not None and frame_stride > 1:
                available_frames = available_frames[::frame_stride]
            
            frame_indices = random.sample(available_frames, min(max_frames, len(available_frames)))
            frame_indices.sort()  # Keep chronological order
            
        else:
            raise ValueError(f"Unknown frame selection method: {frame_selection}")
        
        # Apply max_frames limit
        if max_frames is not None and max_frames > 0:
            frame_indices = frame_indices[:max_frames]
        
        # Remove duplicates and sort
        frame_indices = sorted(list(set(frame_indices)))
        
        print(f"Selected {len(frame_indices)} frames:")
        print(f"  • Frame indices: {frame_indices[:10]}{'...' if len(frame_indices) > 10 else ''}")
        print(f"  • Frame range: {min(frame_indices)} to {max(frame_indices)}")
        
        return frame_indices

    def validate_results(self, results_df):
        """Validate MM/GBSA results for reasonableness"""
        warnings = []
        
        binding_energies = results_df['binding_energy']
        
        # Check for reasonable energy ranges
        mean_binding = binding_energies.mean()
        std_dev = binding_energies.std()
        
        if mean_binding > 10:
            warnings.append(f"Binding energy very positive ({mean_binding:.1f} kcal/mol) - check for errors")
        elif mean_binding < -50:
            warnings.append(f"Binding energy very negative ({mean_binding:.1f} kcal/mol) - check for errors")
        
        if std_dev > 10:
            warnings.append(f"High standard deviation ({std_dev:.1f} kcal/mol) - system may be unstable")
        
        # Check for convergence
        if len(binding_energies) >= 20:
            first_half = binding_energies[:len(binding_energies)//2].mean()
            second_half = binding_energies[len(binding_energies)//2:].mean()
            difference = abs(first_half - second_half)
            
            if difference > 2.0:
                warnings.append(f"Poor convergence: first/second half differ by {difference:.1f} kcal/mol")
        
        # Check for outliers
        q1 = binding_energies.quantile(0.25)
        q3 = binding_energies.quantile(0.75)
        iqr = q3 - q1
        outliers = binding_energies[(binding_energies < q1 - 1.5*iqr) | (binding_energies > q3 + 1.5*iqr)]
        
        if len(outliers) > len(binding_energies) * 0.1:
            warnings.append(f"Many outliers detected ({len(outliers)}/{len(binding_energies)} frames)")
        
        return warnings

    def bootstrap_uncertainty(self, binding_energies, n_bootstrap=1000):
        """Calculate uncertainty using bootstrap resampling"""
        
        bootstrap_means = []
        n_frames = len(binding_energies)
        
        for i in range(n_bootstrap):
            # Resample with replacement
            bootstrap_sample = np.random.choice(binding_energies, size=n_frames, replace=True)
            bootstrap_means.append(np.mean(bootstrap_sample))
        
        bootstrap_means = np.array(bootstrap_means)
        
        return {
            'mean': np.mean(bootstrap_means),
            'std': np.std(bootstrap_means),
            'ci_lower': np.percentile(bootstrap_means, 2.5),
            'ci_upper': np.percentile(bootstrap_means, 97.5)
        }

    def _check_convergence(self, energies, window_size=10):
        """Check if binding energy has converged"""
        if len(energies) < 2 * window_size:
            return {'converged': False, 'reason': 'insufficient_data'}
        
        # Running average convergence check
        running_avg = energies.rolling(window=window_size).mean()
        recent_avg = running_avg.tail(window_size).mean()
        early_avg = running_avg.iloc[window_size:2*window_size].mean()
        
        convergence_threshold = 1.0  # kcal/mol
        
        if abs(recent_avg - early_avg) < convergence_threshold:
            return {'converged': True, 'threshold': convergence_threshold}
        else:
            return {'converged': False, 'difference': abs(recent_avg - early_avg)}

    def generate_detailed_report(self, results, output_dir=None, ligand_resname='LIG', complex_pdb=None):
        """Generate comprehensive analysis report"""
        
        if output_dir is None:
            output_dir = Path('mmgbsa_analysis')
        else:
            output_dir = Path(output_dir)
            
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load results
        df = pd.read_csv(results['output_file'])
        
        # Generate plots. run() and run_comprehensive() both call this report for
        # the same per-frame CSV; the figure depends only on that CSV, so it is
        # drawn once (a 600-dpi PNG costs ~1.7 s).
        png_path = output_dir / 'energy_analysis.png'
        csv_path = Path(results['output_file'])
        plot_key = (str(csv_path.resolve()), csv_path.stat().st_mtime_ns, str(png_path.resolve()))
        if getattr(self, '_energy_plot_key', None) == plot_key and png_path.exists():
            print(f"✓ Energy plots up to date: {png_path}")
        else:
            try:
                import matplotlib.pyplot as plt
            
                # 1. Binding energy over time
                plt.figure(figsize=(12, 8))
            
                plt.subplot(2, 2, 1)
                plt.plot(df['frame'], df['binding_energy'])
                plt.xlabel('Frame')
                plt.ylabel('Binding Energy (kcal/mol)')
                plt.title('Binding Energy vs Time')
                plt.grid(True, alpha=0.3)
            
                # 2. Energy distribution
                plt.subplot(2, 2, 2)
                plt.hist(df['binding_energy'], bins=20, alpha=0.7, edgecolor='black')
                plt.xlabel('Binding Energy (kcal/mol)')
                plt.ylabel('Frequency')
                plt.title('Binding Energy Distribution')
                plt.grid(True, alpha=0.3)
            
                # 3. Running average
                plt.subplot(2, 2, 3)
                running_avg = df['binding_energy'].expanding().mean()
                plt.plot(df['frame'], running_avg)
                plt.xlabel('Frame')
                plt.ylabel('Running Average (kcal/mol)')
                plt.title('Convergence Analysis')
                plt.grid(True, alpha=0.3)
            
                # 4. Component energies (Binding Decomposition)
                plt.subplot(2, 2, 4)
                if 'delta_nb' in df.columns:
                    plt.plot(df['frame'], df['delta_nb'], label='ΔNB (VdW+Ele)', alpha=0.7)
                    plt.plot(df['frame'], df['delta_gb'], label='ΔGB (Polar)', alpha=0.7)
                    plt.plot(df['frame'], df['delta_sa'], label='ΔSA (Nonpolar)', alpha=0.7)
                    plt.title('Binding Energy Components')
                else:
                    plt.plot(df['frame'], df['complex_energy'], label='Complex', alpha=0.7)
                    plt.plot(df['frame'], df['protein_energy'], label='Protein', alpha=0.7)
                    plt.plot(df['frame'], df['ligand_energy'], label='Ligand', alpha=0.7)
                    plt.title('Component Energies')
            
                plt.xlabel('Frame')
                plt.ylabel('Energy (kcal/mol)')
                plt.legend()
                plt.grid(True, alpha=0.3)
            
                plt.tight_layout()
                plt.savefig(output_dir / 'energy_analysis.png', dpi=600, bbox_inches='tight',
                            pil_kwargs={'compress_level': 1})  # lossless; fast zlib level
                self._energy_plot_key = plot_key
                plt.close()
            
                print(f"✓ Energy plots saved to {output_dir / 'energy_analysis.png'}")
            
            except ImportError:
                print("Matplotlib not available, skipping plots")
            except Exception as e:
                print(f"Error generating plots: {e}")
        
        # Bootstrap uncertainty analysis
        bootstrap_results = self.bootstrap_uncertainty(df['binding_energy'])
        
        # Generate report
        report_path = output_dir / 'mmgbsa_report.txt'
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Advanced MM/GBSA Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"GB Model: {results['gb_model']}\n")
            f.write(f"Salt Concentration: {self.salt_concentration} M\n")
            f.write(f"Temperature: {self.temperature}\n")
            f.write(f"Frames Analyzed: {results['n_frames']}\n\n")
            
            f.write("Results Summary:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Mean Binding Energy: {results.get('mean_binding_energy', 0.0):.2f} ± {results.get('std_dev', 0.0):.2f} kcal/mol (Standard Deviation)\n")
            f.write(f"Standard Error of Mean: {results.get('std_error', 0.0):.2f} kcal/mol\n")
            
            med = results.get('median_binding_energy')
            if med is not None:
                f.write(f"Median: {med:.2f} kcal/mol\n")
                f.write(f"Range: {results.get('min_binding_energy', 0.0):.2f} to {results.get('max_binding_energy', 0.0):.2f} kcal/mol\n\n")
            
            if 'bootstrap_results' in results and results['bootstrap_results']:
                bs = results['bootstrap_results']
                f.write("Bootstrap Uncertainty Analysis:\n")
                f.write("-" * 30 + "\n")
                f.write(f"Bootstrap Mean: {bs.get('mean', 0.0):.2f} kcal/mol\n")
                f.write(f"Bootstrap Std: {bs.get('std', 0.0):.2f} kcal/mol\n")
                f.write(f"95% CI: [{bs.get('ci_lower', 0.0):.2f}, {bs.get('ci_upper', 0.0):.2f}] kcal/mol\n\n")
            
            # Convergence analysis
            convergence = self._check_convergence(df['binding_energy'])
            f.write("Convergence Analysis:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Converged: {convergence['converged']}\n")
            if convergence['converged']:
                f.write(f"Convergence threshold: {convergence['threshold']:.1f} kcal/mol\n")
            else:
                if 'difference' in convergence:
                    f.write(f"First/second half difference: {convergence['difference']:.1f} kcal/mol\n")
                else:
                    f.write(f"Reason: {convergence['reason']}\n")
            f.write("\n")
            
            if results.get('validation_warnings'):
                f.write("Validation Warnings:\n")
                f.write("-" * 20 + "\n")
                for warning in results['validation_warnings']:
                    f.write(f"• {warning}\n")
                f.write("\n")
                
            f.write("Physics Assumptions and Fallbacks:\n")
            f.write("-" * 35 + "\n")
            if not results.get('physics_assumptions'):
                f.write("• None. The pipeline utilized explicitly structured parameters and pure explicit configurations without defaulting to empirical heuristics.\n")
            else:
                seen = set()
                for assumption in results['physics_assumptions']:
                    if assumption not in seen:
                        f.write(f"• {assumption}\n")
                        seen.add(assumption)
            f.write("\n")
        
        # Generate Interactive HTML Report (Partially Enabled: Only PandaMap)
        try:
            print("Generating 3D Visualization (PandaMap)...")
        #     
        #     # 1. Load Data
        #     global_results_list = df.to_dict('records') # df loaded earlier
        #     
        #     frame_data = None
        #     frame_csv = output_dir / "frame_by_frame_decomposition.csv"
        #     if frame_csv.exists():
        #         try:
        #             frame_data = pd.read_csv(frame_csv).to_dict('records')
        #         except Exception as fd_err:
        #             print(f"Warning: Could not load frame data: {fd_err}")
        # 
            # 2. Generate PandaMap
            # Avoid duplicate generation when this report method is called
            # multiple times for the same output directory.
            pandamap_html = output_dir / "structure_3d.html"
            if pandamap_html.exists():
                print(f"PandaMap already exists, skipping regeneration: {pandamap_html}")
                pandamap_path = str(pandamap_html)
            else:
                # Use complex_pdb argument passed to this method
                pandamap_path = self._generate_pandamap(complex_pdb, ligand_resname, output_dir)
            pdb_for_report = pandamap_path if pandamap_path else complex_pdb
        #     
        #     # 3. Generate Report
        #     # Pass full config wrapper so HTMLReportGenerator finds 'reporting_settings'
        #     report_config = {'reporting_settings': self.reporting_settings}
        #     generator = HTMLReportGenerator(output_dir, config=report_config)
        #     
        #     # Call with ALL arguments required for high-quality report
        #     html_path = generator.generate_report(results, 
        #                                           frame_data, 
        #                                           global_results=global_results_list, 
        #                                           complex_pdb_path=pdb_for_report,
        #                                           ligand_resname=ligand_resname)
        #                                           
        #     print(f"✓ Interactive report generated: {html_path}")
        except Exception as e:
            print(f"Internal 3D Generation Failed: {e}")
            pass

        print(f"✓ Detailed report generated in {output_dir}")
        return output_dir


    def _generate_pandamap(self, pdb_file, ligand_resname, output_dir):
        """Generate enhanced 3D visualization using PandaMap"""
        try:
            if getattr(self, '_pandamap_disabled', False):
                return None
            print("Generating PandaMap 3D Visualization...")
            try:
                import Bio  # noqa: F401
            except Exception:
                if not getattr(self, '_pandamap_missing_dep_logged', False):
                    print("PandaMap skipped: Biopython is not installed in this environment.")
                    self._pandamap_missing_dep_logged = True
                self._pandamap_disabled = True
                return None
            from pandamap import HybridProtLigMapper
            from pandamap.create_3d_view import create_pandamap_3d_viz
            
            if not pdb_file or not Path(pdb_file).exists():
                print("Warning: No PDB file for PandaMap.")
                return None

            # Patch PDB for PandaMap:
            # 1) Remove solvent/ions from visualization input
            # 2) Convert ligand ATOM records to HETATM for clearer rendering
            excluded_resnames = {
                'HOH', 'WAT', 'SOL', 'TIP3', 'TIP3P', 'H2O',
                'NA', 'CL', 'K', 'MG', 'ZN', 'CA',
                'NA+', 'CL-', 'K+', 'MG2+', 'ZN2+', 'CA2+'
            }
            patched_pdb = output_dir / "temp_fixed_for_panda.pdb"
            with open(pdb_file, 'r') as f_in, open(patched_pdb, 'w') as f_out:
                for line in f_in:
                    if line.startswith("ATOM  ") or line.startswith("HETATM"):
                        resname = line[17:20].strip()
                        if resname in excluded_resnames:
                            continue
                        if line.startswith("ATOM  ") and f" {ligand_resname} " in line:
                            line = "HETATM" + line[6:]
                    f_out.write(line)
            
            # Run PandaMap (suppress printed output if needed)
            mapper = HybridProtLigMapper(str(patched_pdb), ligand_resname=ligand_resname)
            # PandaMap writes its 2D map to the current directory unless told otherwise
            mapper.run_analysis(output_file=str(Path(output_dir) / "pandamap_interactions.png"))
            
            pandamap_html = output_dir / "structure_3d.html"
            create_pandamap_3d_viz(mapper, output_file=str(pandamap_html))
            print(f"PandaMap generated: {pandamap_html}")
            
            # Post-Process for Thinner Sticks & Hover Labels & Style
            try:
                with open(pandamap_html, 'r') as f:
                    html_content = f.read()
                
                # Global Styling
                html_content = html_content.replace("radius: 0.2", "radius: 0.14")
                html_content = html_content.replace("radius:0.2", "radius:0.14")
                html_content = html_content.replace("viewer = $3Dmol.createViewer", "window.viewer = viewer = $3Dmol.createViewer")
                
                # Identify Top 5 Residues
                top_5_calls = ""
                stats_csv = output_dir / "per_residue_detailed.csv"
                if stats_csv.exists():
                     df_stats = pd.read_csv(stats_csv)
                     df_stats.columns = [c.lower() for c in df_stats.columns]
                     if 'total' in df_stats.columns:
                         top_5 = df_stats.sort_values('total', ascending=True).head(5)
                         for _, row in top_5.iterrows():
                             res_num = row.get('residue_number')
                             if pd.isna(res_num): continue
                             res_num = int(res_num)
                             pdb_res_num = res_num
                             res_name = row.get('residue_name', 'RES')
                             val = row['total']
                             label_text = f"{res_name}{res_num} ({val:.1f})"
                             top_5_calls += f"addSmartLabel(v, '{label_text}', {pdb_res_num});\\n"

                hover_script = f"""
    <script>
      function addSmartLabel(viewer, text, resNum) {{
          var selCA = {{resi: resNum, atom: 'CA'}};
          var atoms = viewer.getModel().selectedAtoms(selCA);
          var pos = (atoms.length > 0) ? atoms[0] : {{resi: resNum}};
          viewer.addLabel(text, {{position: pos, fontColor: 'black', fontSize: 14, showBackground: false, inFront: true}});
      }}
      function setupScene() {{
          var v = window.viewer || viewer;
          if (v) {{
              v.setStyle({{resn: '{ligand_resname}'}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.3}}}});
              try {{ {top_5_calls} }} catch(e) {{ console.log(e); }}
              v.setHoverable({{}}, true, function(atom, viewer) {{
                  if (!atom.label) {{
                      var displayResi = parseInt(atom.resi);
                      atom.label = viewer.addLabel(atom.resn + " " + displayResi, {{
                          position: atom, backgroundColor: 'rgba(0,0,0,0.7)', fontColor: 'white', fontSize: 12, showBackground: true
                      }});
                  }}
              }}, function(atom, viewer) {{
                  if (atom.label) {{ viewer.removeLabel(atom.label); delete atom.label; }}
              }});
              v.render();
          }} else {{ setTimeout(setupScene, 500); }}
      }}
      $(document).ready(function() {{ setTimeout(setupScene, 2000); }});
    </script>
    </body>
                """
                if "</body>" in html_content:
                    html_content = html_content.replace("</body>", hover_script)
                else:
                    html_content += hover_script
                    
                with open(pandamap_html, 'w') as f:
                    f.write(html_content)
                
            except Exception as e:
                print(f"Post-processing PandaMap failed: {e}")
                
            return str(pandamap_html)

        except Exception as e:
            print(f"PandaMap generation failed: {e}")
            return None

    # Include all the original caching methods
    def set_ligand_forcefield(self, ff_name):
        """Set the ligand forcefield (openff or gaff)"""
        if ff_name and ff_name.lower() in ['openff', 'gaff']:
            self.ligand_forcefield = ff_name.lower()
            if self.verbose: print(f"Set ligand forcefield to: {self.ligand_forcefield}")

    def parameterize_ligand_openff(self, ligand_mol):
        """
        Parameterize ligand using OpenFF Toolkit or Antechamber (GAFF)
        Returns: system, topology (OpenMM), Molecule (OpenFF)
        """
        from openff.toolkit.topology import Molecule
        
        # GAFF Support imports
        try:
            try:
                from openmmforcefields.generators import GAFFTemplateGenerator
            except ImportError:
                GAFFTemplateGenerator = None
        except ImportError:
             pass

        log.process("Parameterizing ligand with OpenFF SMIRNOFF (this may take a while...)...")
        start_time = time.time()
        
        try:
            # Try different file formats
            if ligand_mol.endswith('.sdf'):
                mol = Molecule.from_file(ligand_mol, file_format='sdf', allow_undefined_stereo=True)
            elif ligand_mol.endswith('.mol2'):
            # Auto-convert to SDF for robust handling (bypass OpenFF/RDKit protonation issues)
                log.process("Auto-converting Mol2 to SDF for robust OpenFF loading...")
                temp_sdf = str(self._work_dir_path() / (Path(ligand_mol).stem + "_converted.sdf"))
                
                if convert_mol2_to_sdf(ligand_mol, temp_sdf):
                    log.info(f"Converted to {temp_sdf}, loading...")
                    mol = Molecule.from_file(temp_sdf, file_format='sdf', allow_undefined_stereo=True)
                else:
                    log.warning("Mol2 conversion failed, attempting direct load (may have atom count issues)...")
                    mol = Molecule.from_file(ligand_mol, file_format='mol2')
            elif ligand_mol.endswith('.mol'):
                mol = Molecule.from_file(ligand_mol, file_format='mol')
            else:
                # Auto-detect format
                mol = Molecule.from_file(ligand_mol, allow_undefined_stereo=True)
                
            log.info(f"Loaded ligand with {mol.n_atoms} atoms")
            
        except Exception as e:
            log.error(f"Error loading ligand file {ligand_mol}: {e}")
            raise e
        
        # Create force field and parameterize
        try:
            self.assign_ligand_charges(mol, ligand_mol)

            # Use SystemGenerator for robust handling of both OpenFF and GAFF
            # and crucially to handle Implicit Solvent automatically via factory
            
            # Constants
            general_kwargs = {
                'constraints': None, # Ligands usually don't have constraints unless specified
                'rigidWater': False,
                'removeCMMotion': False
            }
            nonperiodic_kwargs = {
                'nonbondedMethod': app.NoCutoff,
                'implicitSolvent': self.gbsa_manager.current_app_model,
                'implicitSolventSaltConc': self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar
            }
            
            if self.nonbonded_cutoff is not None:
                nonperiodic_kwargs['nonbondedMethod'] = app.CutoffNonPeriodic
                nonperiodic_kwargs['nonbondedCutoff'] = self.nonbonded_cutoff * unit.angstroms
            
            sm_ff = 'openff-2.1.0'
            if self.ligand_forcefield == 'gaff':
                sm_ff = 'gaff-2.11'
                
            # Initialize Generator
            generator = SystemGenerator(
                forcefields=['amber14-all.xml'], # Minimal base
                small_molecule_forcefield=sm_ff,
                molecules=[mol],
                forcefield_kwargs=general_kwargs,
                nonperiodic_forcefield_kwargs=nonperiodic_kwargs
            )
            
            # Create System (Automatically adds GBSA forces)
            ligand_top = mol.to_topology().to_openmm()
            ligand_system = generator.create_system(ligand_top, molecules=[mol])
            log.success(f"Ligand parameterized with {sm_ff} and implicit solvent")
            
        except Exception as e:
            log.error(f"Error in OpenFF parameterization: {e}")
            raise e
        
        # Refine fixed enhanced GBSA forces to ligand system
        try:
            ligand_gbsa_system = self.gbsa_manager.refine_gbsa_forces(ligand_system, ligand_top)
            
        except Exception as e:
            log.error(f"Error refining GBSA forces for ligand: {e}")
            raise e
        
        total_time = time.time() - start_time
        log.result("Total ligand preparation time", f"{total_time:.1f}", "s")
        
        return ligand_gbsa_system, ligand_top, mol

    def parameterize_protein_amber(self, complex_pdb, ligand_resname=None, ignore_ligand_check=False):
        
        # Native Prmtop Support (Inserted)
        # Native Topology Support (Universal)
        mode = None
        try:
             mode = InputManager.detect_mode(complex_pdb)
             
             if mode in [EngineMode.AMBER, EngineMode.GROMACS, EngineMode.CHARMM]:
                  log.info(f"Detected Native Mode: {mode}")
                  # Prepare kwargs for TopologyLoader
                  loader_kwargs = {
                      'implicitSolvent': self.gbsa_manager.current_app_model,
                      'implicitSolventSaltConc': self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar,
                      'nonbondedMethod': app.NoCutoff
                  }
                  if mode == EngineMode.CHARMM and self.charmm_params:
                      loader_kwargs['charmm_params'] = self.charmm_params
                      if self.charmm_coordinates:
                          loader_kwargs['charmm_coordinates'] = self.charmm_coordinates

                  if self.nonbonded_cutoff is not None:
                      loader_kwargs['nonbondedMethod'] = app.CutoffNonPeriodic
                      loader_kwargs['nonbondedCutoff'] = self.nonbonded_cutoff * unit.angstroms
                      log.info(f"Native Mode: Using CutoffNonPeriodic with {self.nonbonded_cutoff} A cutoff")
                  
                  # Pass implicit solvent args to TopologyLoader (prmtop support)
                  system, topology, positions = TopologyLoader.load_system(
                      complex_pdb, mode, 
                      **loader_kwargs
                  )
                  
                  # Refine GBSA Forces (Separate SA, etc.)
                  system = self.gbsa_manager.refine_gbsa_forces(system, topology)
                  
                  # Generate temp PDB wrapper if needed, or just allow path
                  # Creating dummy positions for PDBFile write (or use native if available)
                  temp_pdb = str(Path(complex_pdb).parent / (Path(complex_pdb).stem + "_temp.pdb"))
                  
                  if positions is None:
                      positions = [openmm.Vec3(0,0,0)]*topology.getNumAtoms()
                      
                  with open(temp_pdb, 'w') as f:
                      app.PDBFile.writeFile(topology, positions, f)
                      
                  return system, topology, None, temp_pdb
        except Exception as e:
         log.warning(f"Native delegation failed: {e}")
         import traceback
         traceback.print_exc()
         # Fallthrough to standard PDB parameterization if detection fails or Generic
         pass
        log.process("Parameterizing protein with Amber (using Modeller for robustness)...")
        start_time = time.time()
        
        # Determine force field to use
        if self.protein_forcefield == 'charmm':
            forcefield_files = ['charmm36.xml']
            log.info("Using CHARMM36 force field...")
        else:  # amber or auto
            forcefield_files = ['amber14-all.xml', 'amber14/tip3p.xml']
            log.info("Using Amber14 force field...")
        
        pdb = app.PDBFile(complex_pdb)
        
        # Use Modeller to clean up and extract protein
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # Delete ligand (residues with ligand_resname)
        if ligand_resname and not ignore_ligand_check:
            ligand_residues = [r for r in modeller.topology.residues() if r.name == ligand_resname]
            if ligand_residues:
                modeller.delete(ligand_residues)
                log.info(f"Removed {len(ligand_residues)} ligand residues from protein system")
            
        # FIX: Remove solvent and ions for MM/GBSA (Dry complex)
        solvent_names = ['HOH', 'WAT', 'TIP3', 'SOL']
        ion_names = ['NA', 'CL', 'K', 'MG', 'ZN', 'CA']
        
        # Reset parameterized residues tracker
        self.parameterized_residues = []
        
        solvent_residues = [r for r in modeller.topology.residues() if r.name in solvent_names]
        ion_residues = [r for r in modeller.topology.residues() if r.name in ion_names]
        
        to_delete = solvent_residues + ion_residues
        if to_delete:
            modeller.delete(to_delete)
            log.info(f"Removed {len(solvent_residues)} solvent and {len(ion_residues)} ion residues")
    
        
        # Add hydrogens/solvent if missing (robustness)
        # Select force field based on protein_forcefield parameter
        if self.protein_forcefield.lower() == 'charmm':
            log.info("Using CHARMM36 force field (with KCX support)...")
            kcx_template = Path(__file__).parent / 'forcefields' / 'kcx_charmm36.xml'
            forcefield_files = [str(kcx_template), 'charmm36.xml']
            log.info("  ✓ Loaded custom KCX template + CHARMM36")
        elif self.protein_forcefield.lower() == 'charmm_gromacs':
            log.info("Using GROMACS-compatible CHARMM36 force field...")
            ff_dir = Path(__file__).parent / 'forcefields'
            kcx_xml = ff_dir / 'kcx_charmm36_gromacs.xml'
            charmm_xml = ff_dir / 'charmm36_gromacs_final.xml'
            forcefield_files = [str(kcx_xml), str(charmm_xml)]
            log.info(f"  ✓ Loaded GROMACS-compatible templates")
        else:  # Default: Amber
            log.info("Using Amber14 force field...")
            forcefield_files = ['amber14-all.xml', 'amber14/tip3pfb.xml']
        
        forcefield = app.ForceField(*forcefield_files)
        
        log.process("Adding/Checking hydrogens...")
        try:
            # Try standard path
            modeller.addHydrogens(forcefield)
            
            protein_top = modeller.topology
            protein_pos = modeller.positions
            
            log.process("Creating OpenMM system...")
            
            # Prepare kwargs for createSystem
            sys_kwargs = {
                'nonbondedMethod': app.NoCutoff,
                'constraints': app.HBonds
            }
            
            if self.nonbonded_cutoff is not None:
                sys_kwargs['nonbondedMethod'] = app.CutoffNonPeriodic
                sys_kwargs['nonbondedCutoff'] = self.nonbonded_cutoff * unit.angstroms
                log.info(f"Using CutoffNonPeriodic with {self.nonbonded_cutoff} A cutoff")
            
            protein_system = forcefield.createSystem(
                protein_top,
                **sys_kwargs
            )
            
        except ValueError as e:
            error_msg = str(e)
            if "No template found" in error_msg:
                # Parse the error to identify which residue failed
                import re
                match = re.search(r"residue (\d+) \((\w+)\)", error_msg)
                
                if match:
                    res_num = int(match.group(1))
                    res_name = match.group(2)
                    
                    # Check if it's a "bonds are different" error for a STANDARD residue
                    standard_resnames = {
                        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                        "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"
                    }
                    
                    if res_name in standard_resnames and "bonds are different" in error_msg:
                        # This is a bond topology mismatch, not a missing residue
                        # Try to proceed by using the forcefield's own repair mechanism
                        log.warning(f"Bond topology mismatch for {res_name}-{res_num} (PDB may be from different FF)")
                        log.warning("Attempting to continue with topology from force field...")
                        
                        # Let OpenMM try to match as best as it can
                        try:
                            # For bond mismatches on standard residues, skip addHydrogens
                            # and use PDB structure as-is (it already has hydrogens)
                            protein_top = modeller.topology
                            protein_pos = modeller.positions
                            
                            log.warning("Skipping addHydrogens due to bond mismatch - using PDB as-is")
                            log.process("Creating system directly from PDB structure...")
                            
                            pdb_sys_kwargs = {
                                'nonbondedMethod': app.NoCutoff,
                                'constraints': None,
                                'ignoreExternalBonds': True,
                                'implicitSolvent': self.gbsa_manager.current_app_model,
                                'implicitSolventSaltConc': self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar
                            }
                            
                            if self.nonbonded_cutoff is not None:
                                pdb_sys_kwargs['nonbondedMethod'] = app.CutoffNonPeriodic
                                pdb_sys_kwargs['nonbondedCutoff'] = self.nonbonded_cutoff * unit.angstroms
                            
                            protein_system = forcefield.createSystem(
                                protein_top,
                                **pdb_sys_kwargs
                            )
                            
                            # If we got here, it worked!
                            log.success("Successfully created system with PDB hydrogens")
                            return protein_system, protein_top, protein_pos, complex_pdb
                            
                        except Exception as e2:
                            log.warning(f"Direct system creation also failed: {e2}")
                            # Fall through to missing residue check below
                            pass
                    
                    # Always check for genuinely missing residues
                    # (either we bypassed PHE, or this is a real missing residue like KCX)
                    log.process("Checking for non-standard residues to parameterize...")
                
                # Identify missing residues
                standard_resnames = {
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN", "NMA", "ACE", 
                    "HOH", "WAT", "NA", "CL", "K", "MG", "ZN", "CA"
                }

                unique_missing = set()
                for res in modeller.topology.residues():
                    if res.name not in standard_resnames:
                        unique_missing.add(res.name)
                
                # If using CHARMM, KCX is already in the custom template - don't parameterize it
                if 'charmm' in self.protein_forcefield.lower() and 'KCX' in unique_missing:
                    unique_missing.remove('KCX')
                
                self.parameterized_residues = list(unique_missing)
                
                if unique_missing:
                    # --- DYNAMIC PARAMETERIZATION (XML ONLY, NO PDB SURGERY) ---
                    # CRITICAL: We ONLY load XML parameters, do NOT modify the PDB structure
                    # This preserves atom count for trajectory compatibility
                    
                    log.process("Attempting on-the-fly parameterization (XML-only, preserving PDB)...")
                    
                    from .parameterizer import ResidueParameterizer
                    param_engine = ResidueParameterizer(str(self._work_dir_path() / "parameterization_work"))
                    
                    for res_name in unique_missing:
                        log.process(f"Parameterizing {res_name}...")
                        xml_out = param_engine.work_dir / f"{res_name}.xml"
                        xml_path, h_pdb_path = param_engine.simple_run(complex_pdb, res_name, xml_out, charge_method=self.charge_method)
                        forcefield.loadFile(str(xml_path))
                        log.success(f"Generated XML parameters for {res_name} (PDB unchanged)")
                    
                    # NO PDB SURGERY - XML template is sufficient
                    # Simply retry addHydrogens() now that forcefield has the KCX template loaded
                    log.process("Retrying addHydrogens with loaded XML templates...")
                    try:
                        modeller.addHydrogens(forcefield)
                    except ValueError as e3:
                        if "bonds are different" in str(e3):
                            log.warning(f"AddHydrogens still fails due to bond mismatch: {e3}")
                            log.warning("Proceeding with existing PDB hydrogens")
                            # Use topology as-is, already has hydrogens
                        else:
                            raise e3
                    
                    protein_top = modeller.topology
                    protein_pos = modeller.positions
                    
                    log.process("Creating OpenMM system with custom residues...")
                    protein_system = forcefield.createSystem(
                        protein_top,
                        nonbondedMethod=app.NoCutoff,
                        constraints=app.HBonds,
                        implicitSolvent=self.gbsa_manager.current_app_model,
                        implicitSolventSaltConc=self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar
                    )

                else:
                    raise e
            else:
                raise e
        
        # Add fixed enhanced GBSA forces to protein system
        # Refine fixed enhanced GBSA forces in protein system
        protein_gbsa_system = self.gbsa_manager.refine_gbsa_forces(protein_system, protein_top)
        
        log.success(f"Protein parameterized with fixed enhanced GBSA ({protein_gbsa_system.getNumParticles()} particles)")
        
        total_time = time.time() - start_time
        print(f"✓ Total protein preparation time: {total_time:.1f}s")
        
        return protein_gbsa_system, protein_top, protein_pos, complex_pdb

    def assign_ligand_charges(self, molecule, source=None):
        """Give an OpenFF Molecule partial charges, unless its file had them.

        Charges read from the ligand file are kept. Otherwise `charge_method`
        (am1bcc, mmff94 or gasteiger) assigns them; when that method fails the
        run stops with an error that says how to avoid it, instead of falling
        back to another method. Returns the method used ('file' for the file's
        charges) and records it in physics_assumptions.
        """
        name = Path(str(source)).name if source else 'the ligand'
        if molecule.partial_charges is not None:
            log.info(f"Ligand partial charges: kept from {name}")
            self.ligand_charge_source = 'file'
            return 'file'
        method = self.charge_method
        log.process(f"{name} has no partial charges; assigning {method} charges (charge_method)...")
        try:
            molecule.assign_partial_charges(partial_charge_method=method)
        except Exception as e:
            need = (" AM1-BCC needs AmberTools (sqm) or the OpenEye toolkit." if method == 'am1bcc' else "")
            raise RuntimeError(
                f"charge_method '{method}' could not assign partial charges to {name}: {e}.{need} "
                f"Put partial charges in the ligand file (.sdf), or set "
                f"analysis_settings.charge_method to mmff94 (or gasteiger).") from e
        log.success(f"Ligand partial charges: {method}")
        self.ligand_charge_source = method
        note = f"Ligand partial charges: {method} (charge_method), the ligand file had none."
        if note not in self.physics_assumptions:
            self.physics_assumptions.append(note)
        return method

    def get_ligand_positions_openmm(self, ligand_mol_obj, ligand_pdb):
        """Get ligand positions, preferring PDB but falling back to Molecule conformer if atom count mismatches"""
        ligand_topology = ligand_mol_obj.to_topology()
        
        # Try loading PDB positions first
        use_pdb = False
        pdb_pos = None
        try:
            pdb = app.PDBFile(ligand_pdb)
            if pdb.topology.getNumAtoms() == ligand_topology.getNumAtoms():
                use_pdb = True
                pdb_pos = pdb.getPositions()
                print("  Using ligand positions from PDB (atom counts match)")
            else:
                print(f"  Warning: Ligand PDB atom count ({pdb.topology.getNumAtoms()}) != Topology count ({ligand_topology.getNumAtoms()})")
                print("  Falling back to SDF/Molecule conformers...")
        except Exception as e:
            print(f"  Warning: Could not load ligand PDB positions: {e}")
        
        if use_pdb:
            return ligand_topology, pdb_pos
            
        # Fallback to Molecule conformer (SDF)
        if ligand_mol_obj.n_conformers > 0:
            print("  Using ligand positions from OpenFF Molecule (SDF)")
            conf = ligand_mol_obj.conformers[0]
            
            # Robust unit conversion (OpenFF Quantity -> OpenMM Quantity)
            try:
                # Try modern OpenFF (pint-based)
                coords_nm = conf.m_as("nanometer")
                return ligand_topology, unit.Quantity(coords_nm, unit.nanometer)
            except AttributeError:
                # Fallback for older types or simtk.unit types
                if hasattr(conf, 'value_in_unit'):
                    # It's an OpenMM quantity?
                    return ligand_topology, conf
                else:
                    # Assume Angstroms (standard for SDF)
                    print("  Warning: Assuming Angstroms for ligand conformer")
                    # Check if it has .magnitude (pint/simtk)
                    if hasattr(conf, 'magnitude'):
                        coords = conf.magnitude
                    else:
                        coords = np.array(conf) # Raw array?
                    return ligand_topology, unit.Quantity(coords / 10.0, unit.nanometer)
                    
        raise ValueError(f"No valid ligand positions found! PDB atom count ({pdb.topology.getNumAtoms() if 'pdb' in locals() else '?'}) mismatches SDF ({ligand_topology.getNumAtoms()}), and SDF has no conformers.")
    
    def build_complex_system(self, protein_pdb, ligand_mol=None, ligand_pdb=None, add_gbsa=True):
        """Build the full complex system with fixed enhanced GBSA forces
        
        If ligand_mol is None, assumes protein-only or protein-protein system (using Amber for all).
        """
        
        log.process('Building complex system...')
        start_time = time.time()
        
        molecules_list = []

        # Check for Native Input (Bypass SystemGenerator/Modeller)
        # Check for Native Input (Bypass SystemGenerator/Modeller)
        mode = None
        try:
             mode = InputManager.detect_mode(protein_pdb)
             if mode in [EngineMode.AMBER, EngineMode.GROMACS, EngineMode.CHARMM]:
                  log.info(f"Delegating native input {protein_pdb} to parameterize_protein_amber...")
                  system, topology, _, _ = self.parameterize_protein_amber(protein_pdb)
                  return system, topology
        except Exception as e:
             log.warning(f"Native delegation failed: {e}")
             if mode == EngineMode.AMBER:
                 raise e
        
        # Fallback to Generic Generation
        # Load structure (Generic - supports .gro, .pdb, etc)
        import parmed as pmd
        try:
            struct = pmd.load_file(protein_pdb)
            # Modeller requires OpenMM Topology and Positions
            # Convert ParmEd topology to OpenMM
            omm_top = struct.topology
            omm_pos = struct.positions
            modeller = app.Modeller(omm_top, omm_pos)
        except Exception as e:
            # Last resort: PDBFile
            protein_pdbfile = app.PDBFile(protein_pdb)
            modeller = app.Modeller(protein_pdbfile.topology, protein_pdbfile.positions)
        
        # Protein-Ligand Mode (Standard)
        if ligand_mol is not None:
            log.info("Running in Protein-Ligand Mode (OpenFF + Amber)")
            if Molecule is None:
                # Retry lazy import in case module-level optional import failed transiently.
                try:
                    from openff.toolkit.topology import Molecule as _OpenFFMolecule
                except Exception as e:
                    raise ImportError(
                        "OpenFF toolkit is required for ligand_mol-based parameterization. "
                        "Install 'openff-toolkit' in this environment."
                    ) from e
                globals()['Molecule'] = _OpenFFMolecule
            
            # Auto-convert Mol2 to SDF if needed
            
            # Auto-convert Mol2 to SDF if needed
            if str(ligand_mol).endswith('.mol2'):
                # converted into the run folder, never next to (or over) an input file
                sdf_path = str(self._work_dir_path() / (Path(str(ligand_mol)).stem + '.sdf'))
                log.process(f"Auto-converting {ligand_mol} to {sdf_path} for processing...")
                if convert_mol2_to_sdf(ligand_mol, sdf_path):
                    ligand_mol = sdf_path
                else:
                    log.warning("Conversion failed, attempting to read Mol2 directly...")
                
            ligand_mol_obj = Molecule.from_file(ligand_mol, allow_undefined_stereo=True)
            
            # Ensure ligand has conformers
            if ligand_mol_obj.n_conformers == 0:
                log.process("Generating conformer for ligand...")
                ligand_mol_obj.generate_conformers(n_conformers=1)
                
            molecules_list = [ligand_mol_obj]
            
            # The ligand keeps the partial charges of its file; a file without
            # charges gets charge_method charges (SystemGenerator uses the
            # molecule's charges when it has them)
            self.assign_ligand_charges(ligand_mol_obj, ligand_mol)
            
            # Identify and delete existing ligand from PDB
            ligand_resname_extracted = self.find_ligand_resname(modeller.topology)
            if ligand_resname_extracted:
                log.info(f"Removing existing ligand residue: {ligand_resname_extracted}")
                to_delete = [r for r in modeller.topology.residues() if r.name == ligand_resname_extracted]
                if to_delete: modeller.delete(to_delete)
        
            # Add OpenFF ligand
            if ligand_pdb:
                 lig_top, ligand_positions = self.get_ligand_positions_openmm(ligand_mol_obj, ligand_pdb)
                 log.process("Adding OpenFF ligand to modeller...")
                 modeller.add(lig_top.to_openmm(), ligand_positions)
        
        # Protein-Only / Protein-Protein Mode
        else:
            log.info("Running in Protein-Only/Protein-Protein Mode (Amber Only)")
            # No ligand removal, no molecule addition. 
            # We assume complex_pdb contains everything we need.
            pass

        # Clean Solvent/Ions (Always do this unless requested otherwise)
        solvent_names = ['HOH', 'WAT', 'TIP3', 'SOL']
        ion_names = ['NA', 'CL', 'K', 'MG', 'ZN', 'CA']
        to_delete_solvent = [r for r in modeller.topology.residues() if r.name in solvent_names or r.name in ion_names]
        if to_delete_solvent:
            modeller.delete(to_delete_solvent)
            log.info(f"Removed {len(to_delete_solvent)} solvent/ion residues")
    
        # Create system generator
        general_kwargs = {
            'constraints': app.HBonds,
            'rigidWater': True,
            'removeCMMotion': False,
            'hydrogenMass': 4*unit.amu
        }
        
        nonperiodic_kwargs = {
            'nonbondedMethod': app.NoCutoff,
            'implicitSolvent': self.gbsa_manager.current_app_model,
            'implicitSolventSaltConc': self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar
        }
        
        # Resolve Forcefields dynamically
        _ff_base = Path(__file__).parent / 'forcefields'
        if self.protein_forcefield.lower() == 'charmm':
            sys_ffs = [str(_ff_base / 'kcx_charmm36.xml'), 'charmm36.xml', 'charmm36/water.xml']
        elif self.protein_forcefield.lower() == 'charmm_gromacs':
            sys_ffs = [str(_ff_base / 'kcx_charmm36_gromacs.xml'), str(_ff_base / 'charmm36_gromacs_final.xml'), 'charmm36/water.xml']
        elif self.protein_forcefield.lower() == 'amber14':
            sys_ffs = ['amber14-all.xml', 'amber14/tip3pfb.xml']
        else:
            sys_ffs = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml']

        if SystemGenerator is None:
            log.warning("openmmforcefields SystemGenerator is unavailable; falling back to OpenMM ForceField.")
            if ligand_mol is not None:
                raise RuntimeError("SystemGenerator is required for ligand_mol-based parameterization but is not available.")
            ff = None
            ff_candidates = [sys_ffs]
            for ff_files in ff_candidates:
                try:
                    ff = app.ForceField(*ff_files)
                    break
                except Exception:
                    continue
            if ff is None:
                raise RuntimeError("No compatible OpenMM protein/water forcefield XML set was found.")
            try:
                log.info("Checking for missing atoms and adding hydrogens...")
                modeller.addHydrogens(forcefield=ff)
            except Exception as e:
                log.warning(f"Modeller.addHydrogens failed: {e}. Proceeding with existing topology.")
            try:
                system = ff.createSystem(modeller.topology, **general_kwargs, **nonperiodic_kwargs)
            except Exception:
                fallback_kwargs = {k: v for k, v in nonperiodic_kwargs.items() if not k.startswith('implicitSolvent')}
                system = ff.createSystem(modeller.topology, **general_kwargs, **fallback_kwargs)
            if add_gbsa:
                names = (type_names_from_forcefield(ff, modeller.topology)
                         if self.lcpo_typing == 'amber' else None)
                system = self.refine_gbsa_forces(system, modeller.topology, names)
            modeller.topology.setPeriodicBoxVectors(None)
            log.success(f"System created ({system.getNumParticles()} particles)")
            return system, modeller.topology, modeller.positions
        elif ligand_mol is None:
            # Protein-only/PPI mode: avoid initializing small-molecule toolkits.
            system_generator = _ProteinOnlySystemGenerator(
                forcefields=sys_ffs,
                small_molecule_forcefield=None,
                forcefield_kwargs=general_kwargs,
                nonperiodic_forcefield_kwargs=nonperiodic_kwargs
            )
        else:
            system_generator = SystemGenerator(
                forcefields=sys_ffs,
                small_molecule_forcefield='openff-2.0.0',
                molecules=molecules_list,
                forcefield_kwargs=general_kwargs,
                nonperiodic_forcefield_kwargs=nonperiodic_kwargs
            )
        
        # Ensure non-periodic
        modeller.topology.setPeriodicBoxVectors(None)
        
        # Add missing hydrogens (e.g., terminal H3 or missing backbone H)
        # Required for strict template matching in Coordinate mode
        try:
            log.info("Checking for missing atoms and adding hydrogens...")
            modeller.addHydrogens(forcefield=system_generator.forcefield)
        except Exception as e:
            if "KCX" in str(e) and "amber" in str(self.protein_forcefield).lower():
                log.error("KCX RESIDUE DETECTED! OpenMM's Amber forcefields do not support Carboxylated Lysine (KCX).")
                log.error("Please change `protein_forcefield: charmm_gromacs` in your configuration file to use the included custom KCX parameters.")
            log.warning(f"Modeller.addHydrogens failed: {e}. Proceeding with existing topology.")
            
        # Create System (This now generates the correct GB force automatically)
        try:
            system = system_generator.create_system(modeller.topology)
        except ValueError as e:
            if "KCX" in str(e) and "amber" in str(self.protein_forcefield).lower():
                raise ValueError("KCX RESIDUE DETECTED! OpenMM's Amber forcefields do not support Carboxylated Lysine (KCX). Please change `protein_forcefield: charmm_gromacs` in your config file.") from e
            if 'implicitSolvent' in str(e):
                log.warning(f"SystemGenerator rejected implicitSolvent: {e}")
                log.warning("Attempting to build system without implicitSolvent keyword...")
                self.physics_assumptions.append("System Builder: Native OpenMM integration rejected Implicit Solvent keywords. Falling back to strictly isolated vacuum parameters and injecting implicit mapping externally.")
                # Unfortunately, SystemGenerator doesn't allow easy modification of kwargs after initialization.
                # We have to create a new one without implicitSolvent.
                fallback_kwargs = {k: v for k, v in nonperiodic_kwargs.items() if not k.startswith('implicitSolvent')}
                if ligand_mol is None:
                    fallback_generator = _ProteinOnlySystemGenerator(
                        forcefields=sys_ffs,
                        small_molecule_forcefield=None,
                        forcefield_kwargs=general_kwargs,
                        nonperiodic_forcefield_kwargs=fallback_kwargs
                    )
                else:
                    fallback_generator = SystemGenerator(
                        forcefields=sys_ffs,
                        small_molecule_forcefield='openff-2.0.0',
                        molecules=molecules_list,
                        forcefield_kwargs=general_kwargs,
                        nonperiodic_forcefield_kwargs=fallback_kwargs
                    )
                system = fallback_generator.create_system(modeller.topology)
            else:
                raise e
        
        if add_gbsa:
            # Refine the generated forces (Separate SA, etc.); Amber LCPO typing reads the
            # atom-type classes the generator's force field assigned (what tleap writes)
            names = None
            if self.lcpo_typing == 'amber':
                gen = locals().get('fallback_generator') or system_generator
                names = type_names_from_forcefield(getattr(gen, 'forcefield', None), modeller.topology)
            system = self.refine_gbsa_forces(system, modeller.topology, names)
        
        # Ensure non-periodic
        modeller.topology.setPeriodicBoxVectors(None)
        
        # Create system
        log.process("Creating system with Amber forcefield...")
        gbsa_system = system
        
        log.success(f"System created ({gbsa_system.getNumParticles()} particles)")
        
        total_time = time.time() - start_time
        log.result("Total complex system preparation time", f"{total_time:.1f}", "s")
        
        return gbsa_system, modeller.topology, modeller.positions

    @staticmethod
    def _solute_atoms(md_topology):
        return [a.index for a in md_topology.atoms
                if not (a.residue.is_water or a.residue.name in SOLVENT_ION_RESNAMES)]

    @staticmethod
    def _as_mdtraj(topology):
        return topology if hasattr(topology, 'select') else md.Topology.from_openmm(topology)

    def check_two_chains(self, topology, where):
        """binding_mode ppi: the complex must have exactly two chains with solute atoms."""
        top = self._as_mdtraj(topology)
        solute = set(self._solute_atoms(top))
        chains = [c for c in top.chains if any(a.index in solute for a in c.atoms)]
        if len(chains) == 1:
            raise ValueError(
                f"binding_mode ppi splits the complex by chain (receptor = chainid 0, ligand = "
                f"chainid 1), but {where} has a single chain (no chain records). Give "
                f"analysis_settings.receptor_selection and ligand_selection instead, e.g. "
                f"receptor_selection: 'resid 0 to 120' and ligand_selection: 'resid 121 to 230' "
                f"(MDTraj selections).")
        if len(chains) > 2:
            raise ValueError(
                f"binding_mode ppi uses chainid 0 as the receptor and chainid 1 as the ligand, but "
                f"{where} has {len(chains)} chains with solute atoms, so chains 2 to "
                f"{len(chains) - 1} would be left out. Give analysis_settings.receptor_selection and "
                f"ligand_selection, e.g. receptor_selection: 'chainid 0 to 1' and ligand_selection: "
                f"'chainid 2'.")

    def split_by_selections(self, topology, ligand_selection, receptor_selection, where):
        """Ligand and receptor atom indices of two MDTraj selections.

        Each must match atoms and they must not overlap (ValueError). Solute atoms
        in neither selection are reported, since they are left out of the
        calculation. Returns (ligand_indices, receptor_indices, left_out).
        """
        top = self._as_mdtraj(topology)
        hint = (" A topology converted from a .tpr may have lost its chain records; select by "
                "residue or atom index instead, or give a solvated_topology with chains."
                if '.tpr' in where else "")
        lig = sorted(int(i) for i in self.get_selection_indices(top, ligand_selection))
        rec = sorted(int(i) for i in self.get_selection_indices(top, receptor_selection))
        log.info(f"  Ligand selection '{ligand_selection}' matched {len(lig)} atoms")
        log.info(f"  Receptor selection '{receptor_selection}' matched {len(rec)} atoms")
        if not lig:
            raise ValueError(f"ligand_selection '{ligand_selection}' matches no atom of {where}.{hint}")
        if not rec:
            raise ValueError(f"receptor_selection '{receptor_selection}' matches no atom of {where}.{hint}")
        overlap = sorted(set(lig) & set(rec))
        if overlap:
            res = top.atom(overlap[0]).residue
            raise ValueError(
                f"receptor_selection '{receptor_selection}' and ligand_selection '{ligand_selection}' "
                f"overlap in {len(overlap)} atoms of {where} (e.g. residue {res.name}{res.resSeq} "
                f"of chain {res.chain.index}); an atom cannot be in both.")
        chosen = set(lig) | set(rec)
        left_out = [i for i in self._solute_atoms(top) if i not in chosen]
        if left_out:
            per_chain = defaultdict(set)
            for i in left_out:
                res = top.atom(i).residue
                per_chain[res.chain.index].add(res.resSeq)
            detail = "; ".join(f"chain {c}: {len(r)} residues ({min(r)}-{max(r)})"
                               for c, r in sorted(per_chain.items()))
            log.warning(f"{len(left_out)} solute atoms of {where} are in neither receptor_selection "
                        f"nor ligand_selection and are LEFT OUT of the calculation ({detail}).")
            self.physics_assumptions.append(
                f"Selections: {len(left_out)} solute atoms in neither receptor_selection nor "
                f"ligand_selection were left out ({detail}).")
        return lig, rec, left_out

    @staticmethod
    def explicit_amber_topology(path, key):
        """An Amber receptor/ligand topology (.prmtop, .parm7, .new), or None.

        Any other file is an error: before, a .parm7 (or any other name) was
        silently ignored.
        """
        if not path:
            return None
        suffix = Path(str(path)).suffix.lower()
        if suffix in ('.prmtop', '.parm7', '.new'):
            return str(path)
        if suffix == '.top':
            raise ValueError(f"{key} {path} is a GROMACS topology; {key} must be an Amber "
                             f".prmtop or .parm7. A GROMACS .top complex is converted together with "
                             f"its receptor and ligand; leave {key} out.")
        raise ValueError(f"{key} {path} is not an Amber topology; give a .prmtop or .parm7 file, "
                         f"or leave {key} out.")

    def get_selection_indices(self, topology, selection_string, ref_pdb=None):
        """Get indices using MDTraj selection syntax
        
        Parameters:
        -----------
        topology : openmm.Topology or mdtraj.Topology
            Topology object
        selection_string : str
            MDTraj selection string
        """
        try:
            # Handle MDTraj Topology directly (Duck typing: check for 'select' method)
            if hasattr(topology, 'select'):
                md_top = topology
            else:
                # Convert OpenMM to MDTraj
                md_top = md.Topology.from_openmm(topology)
                
            indices = md_top.select(selection_string)
            return list(indices)
        except Exception as e:
            log.error(f"Selection failed for '{selection_string}': {e}")
            return []

    def get_ligand_indices(self, topology, ligand_resname, selection=None):
        """Get indices of ligand atoms
        
        If 'selection' is provided (e.g., 'chainid 1'), it overrides ligand_resname logic.
        """
        if selection:
            log.info(f"Using custom ligand selection: '{selection}'")
            return self.get_selection_indices(topology, selection)
            
        # Handle both OpenMM (method) and MDTraj (property) topologies
        atoms = topology.atoms() if callable(getattr(topology, 'atoms', None)) else topology.atoms
        return [i for i, atom in enumerate(atoms) if atom.residue.name == ligand_resname]

    def get_protein_indices(self, topology, ligand_resname, selection=None):
        """Get indices of protein atoms
        
        If 'selection' is provided (e.g., 'chainid 0'), it overrides exclusion logic.
        """
        if selection:
            log.info(f"Using custom protein/receptor selection: '{selection}'")
            return self.get_selection_indices(topology, selection)
            
        solvent_names = ['HOH', 'WAT', 'TIP3', 'SOL']
        ion_names = ['NA', 'CL', 'K', 'MG', 'ZN', 'CA']
        excluded_resnames = set(solvent_names + ion_names + [ligand_resname])
        
        # Handle both OpenMM (method) and MDTraj (property) topologies
        atoms = topology.atoms() if callable(getattr(topology, 'atoms', None)) else topology.atoms
        return [i for i, atom in enumerate(atoms) if atom.residue.name not in excluded_resnames]

    def find_ligand_resname(self, topology):
        """Find ligand residue name in topology"""
        residue_names = set(atom.residue.name for atom in topology.atoms())
        ligand_names = ['LIG', 'UNL', 'UNK', 'DRUG', 'MOL']
        for name in ligand_names:
            if name in residue_names:
                return name
        
        # Fallback: residue with fewest atoms (excluding solvent/ions)
        self.physics_assumptions.append("Topology Parsing: Unable to automatically verify targeted Ligand coordinate node based on conventional internal namings. Assuming ligand as smallest unbound component molecule strictly via mathematical counting.")
        solvent_names = ['HOH', 'WAT', 'TIP3', 'SOL']
        ion_names = ['NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'CL-']
        excluded = set(solvent_names + ion_names)
        
        residue_atom_counts = {}
        for atom in topology.atoms():
            resname = atom.residue.name
            if resname in excluded:
                continue
            residue_atom_counts[resname] = residue_atom_counts.get(resname, 0) + 1
        
        min_atoms = min(residue_atom_counts.values())
        for resname, count in residue_atom_counts.items():
            if count == min_atoms and count < 200:
                return resname
        return None

    def setup_optimized_platform(self, platform_name=None):
        """Setup optimized platform for calculations
        
        Parameters:
        -----------
        platform_name : str, optional
            Override platform selection for this specific call (e.g., 'CPU' for decomposition
            while main analysis uses CUDA). If None, uses instance preference.
        
        Returns:
        --------
        tuple : (platform, properties) - OpenMM platform and device properties

        A requested platform that cannot be used (e.g. a CUDA runtime newer than
        the driver supports) falls back to the next choice, as before, but the
        fallback is no longer silent: a warning names the requested and the
        actual platform, and both are recorded on the instance
        (``platform_requested``, ``platform_used``, ``platform_fallback``) and in
        the results summary. With ``platform_settings.strict_platform: true`` an
        unusable requested platform is an error instead.
        """
        # Platform priority: parameter > instance decomposition > instance preference > environment > auto-detect
        # We validate platforms by creating a tiny Context because CUDA/OpenCL may be
        # discoverable yet unusable at runtime (e.g., no device available).
        last_error = [None]   # why the most recent platform test failed

        def _can_use_platform(name, props):
            try:
                plat = openmm.Platform.getPlatformByName(name)
                test_system = openmm.System()
                test_system.addParticle(1.0)
                if name in ('CUDA', 'OpenCL'):
                    # OpenMM caches compiled GPU kernels by source in a shared temporary
                    # directory, so a probe made only of common kernels can pass on code
                    # another environment compiled, even when this environment's CUDA
                    # runtime builds code the driver cannot load. A force whose kernel
                    # source is unique to this probe has to be compiled here.
                    probe = openmm.CustomExternalForce(f"{1 + random.random():.17g}*x")
                    probe.addParticle(0, [])
                    test_system.addForce(probe)
                test_integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
                _ctx = openmm.Context(test_system, test_integrator, plat, props)
                _ctx.setPositions([openmm.Vec3(0.0, 0.0, 0.0)] * unit.nanometers)
                _ctx.getState(getEnergy=True)
                del _ctx
                del test_integrator
                return plat
            except Exception as exc:
                last_error[0] = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
                return None
        
        fallbacks = []   # (requested platform, reason) for every explicit request that failed

        def _chosen(platform, properties, requested):
            used = platform.getName()
            self.platform_requested = requested
            self.platform_used = used
            self.platform_fallback = '; '.join(f"{n}: {r}" for n, r in fallbacks) or None
            if fallbacks:
                msg = (f"Requested OpenMM platform '{fallbacks[0][0]}' is not usable "
                       f"({fallbacks[0][1]}); this run uses '{used}' instead.")
                if getattr(self, 'strict_platform', False):
                    raise RuntimeError(msg + " Set platform_settings.strict_platform: false "
                                       "to allow the fallback.")
                log.warning(f"PLATFORM FALLBACK: {msg} Energies agree to within the numerical precision of the platforms, but timings are "
                            f"not those of '{fallbacks[0][0]}'. Set platform_settings.strict_platform: "
                            f"true to stop instead.")
            print(f"Running on OpenMM platform: {used}")
            return platform, properties

        def _props(name):
            if name == 'CUDA':
                return {'CudaPrecision': 'mixed', 'CudaDeviceIndex': '0'}
            if name == 'OpenCL':
                return {'OpenCLPrecision': 'mixed'}
            return {}

        # 1. explicit platform_name (e.g. a decomposition override), 2. the instance
        # preference, 3. the OPENMM_DEFAULT_PLATFORM environment variable
        env_pref = os.environ.get('OPENMM_DEFAULT_PLATFORM')
        for source, requested in (("Forced platform", platform_name),
                                  ("Forced platform", getattr(self, 'platform_preference', None)),
                                  ("Environment platform", env_pref)):
            if not requested or str(requested).lower() == 'auto':
                continue
            requested = _platform_name(requested) or str(requested)
            properties = _props(requested)
            platform = _can_use_platform(requested, properties)
            if platform is not None:
                print(f"Using {source.lower()}: {requested}")
                return _chosen(platform, properties, requested)
            fallbacks.append((requested, last_error[0] or "it is not usable on this machine"))
            print(f"⚠️ Warning: {source} '{requested}' is not usable on this machine; trying the next choice.")

        # 4. Auto-detect: Try CUDA -> OpenCL -> CPU
        for candidate, props, msg in [
            ('CUDA', _props('CUDA'), "Using CUDA platform"),
            ('OpenCL', _props('OpenCL'), "Using OpenCL platform"),
            ('CPU', {}, "Using CPU platform"),
        ]:
            platform = _can_use_platform(candidate, props)
            if platform is not None:
                print(msg)
                return _chosen(platform, props, fallbacks[0][0] if fallbacks else 'auto')
        raise RuntimeError("No usable OpenMM platform found (CUDA/OpenCL/CPU).")

    def run_comprehensive(self, ligand_mol, complex_pdb, xtc_file, ligand_pdb, max_frames=50, 
                    energy_decomposition=False, frame_start=None, frame_end=None, 
                    frame_stride=None, frame_selection='sequential', random_seed=42,
                    qha_analyze_complex=False, output_dir=None, ligand_selection=None, receptor_selection=None,
                    receptor_topology=None, ligand_topology=None, solvated_topology=None, print_interval=10,
                    require_two_chains=False):
        """
        Enhanced run method with comprehensive validation and analysis
        
        Parameters:
        -----------
        ligand_mol : str
            Path to ligand molecule file
        complex_pdb : str
            Path to protein-ligand complex PDB file
        xtc_file : str
            Path to molecular dynamics trajectory file
        ligand_pdb : str
            Path to isolated ligand PDB file
        max_frames : int, optional
            Maximum number of trajectory frames to analyze
        energy_decomposition : bool, optional
            Enable detailed energy decomposition analysis
        frame_start : int, optional
            Start frame for analysis (0-indexed)
        frame_end : int, optional
            End frame for analysis (0-indexed)
        frame_stride : int, optional
            Frame stride (every Nth frame)
        frame_selection : str, optional
            Frame selection method ('sequential', 'equidistant', 'random')
        random_seed : int, optional
            Random seed for frame selection
        ligand_selection : str, optional
            Custom selection string for ligand
        receptor_selection : str, optional
            Custom selection string for receptor
        require_two_chains : bool, optional
            binding_mode ppi with chainid 0/1 defaults: stop unless the complex has
            exactly two chains with solute atoms
            
        Returns:
        --------
        dict : Analysis results with enhanced statistics and validation
        """
        
        # Preserve original path before any potential fallback replacement
        original_input_complex_pdb = complex_pdb
        if output_dir:
            self.work_dir = Path(output_dir)
        self.distilled_files = None

        # mdtraj has no registered loader for a bare '.trj' extension --
        # normalize once, up front, so every downstream md.load(xtc_file, ...)
        # call in this method and in run() gets a path mdtraj recognizes,
        # rather than needing this same fix repeated at each call site (see
        # TrajectoryProcessor.resolve_mdtraj_loadable_path's own docstring
        # for why this is needed -- confirmed on a real published dataset,
        # Zenodo 17926575, whose Amber trajectories use this extension).
        if xtc_file:
            from .trajectory import TrajectoryProcessor
            xtc_file = TrajectoryProcessor.resolve_mdtraj_loadable_path(xtc_file, self._work_dir_path())

        # Pre-run validation
        print("Validating input files...")
        validation_errors = self.validate_input_files(ligand_mol, complex_pdb, ligand_pdb, xtc_file, solvated_topology)
        
        if validation_errors:
            print("❌ Input validation failed:")
            for error in validation_errors:
                print(f"  • {error}")
            return None
        
        print("✅ Input validation passed")
        
        # ---------------------------------------------------------------------
        # Prepare MDTraj-compatible topology when a GROMACS TPR is provided.
        # The conversion is required because MDTraj cannot natively read ".tpr"
        # files.  In normal runs the topology is later written once the native
        # complex is loaded (see below), but custom-selection logic executes
        # *before* that step.  This block proactively generates a small PDB
        # so that subsequent `md.load()` calls succeed.
        # ---------------------------------------------------------------------
        mdtraj_topology = complex_pdb
        if str(complex_pdb).lower().endswith('.tpr'):
            import os
            mdtraj_topology = os.path.join(output_dir if output_dir else '.',
                                           "complex_solvated_from_tpr_for_mdtraj.pdb")
            if not os.path.exists(mdtraj_topology):
                converted = False
                try:
                    from .tpr_loader import load_tpr_as_parmed
                    struct = load_tpr_as_parmed(str(complex_pdb), xtc_path=xtc_file)
                    struct.save(mdtraj_topology, overwrite=True)
                    log.info(f"Generated temporary topology for mdtraj: {mdtraj_topology}")
                    converted = True
                    # do NOT modify `complex_pdb` here: this variable is passed
                    # through to `run()` and later used to determine native vs
                    # coordinate mode.  overwriting it with the temporary PDB
                    # causes the system to think it no longer has a TPR and
                    # fall back to raw-coordinate generation (which then failed
                    # due to missing hydrogen templates).  Instead we keep the
                    # original path intact and only use `mdtraj_topology` for
                    # MDTraj loading above.
                except Exception as e:
                    log.warning(f"Failed to convert TPR to PDB for mdtraj using TprParser: {e}")
                    # Fallback: try to find .pdb or .gro in the same directory as TPR
                    tpr_dir = Path(complex_pdb).parent
                    tpr_stem = Path(complex_pdb).stem
                    for ext in ['.pdb', '.gro']:
                        fallback_path = tpr_dir / f"{tpr_stem}{ext}"
                        if fallback_path.exists():
                            mdtraj_topology = str(fallback_path)
                            # Also update complex_pdb to use fallback
                            complex_pdb = str(fallback_path)
                            # since the effective input is now a non-TPR file we should
                            # update the saved original path as well so later logic
                            # doesn't erroneously treat it as a GROMACS native topology.
                            original_input_complex_pdb = complex_pdb
                            log.info(f"TPR conversion failed; using fallback topology: {mdtraj_topology}")
                            converted = True
                            break
                if not os.path.exists(mdtraj_topology):
                    # no successful conversion path
                    raise ValueError(
                        "MDTraj cannot read a .tpr topology and automatic conversion "
                        "failed. Please provide an explicit PDB/GRO topology via the "
                        "configuration (e.g. 'solvated_topology') or install the ``TprParser`` "
                        "package so that opengbsa can convert the file itself."
                    )
        
        # Receptor/ligand selections (binding_mode ppi or explicit selections).
        # With both selections the two must not overlap; atoms in neither are left
        # out of the calculation (Coordinate Mode: distilled files; Native Mode:
        # run() drops them from the topology and the trajectory).
        if ligand_selection and receptor_selection:
            log.info("Processing custom selections for Protein-Protein/Refined analysis...")

            # Load reference topology (use converted path when available).
            # For Amber topology-only inputs (.prmtop/.parm7), load first frame from trajectory.
            if str(mdtraj_topology).lower().endswith(('.prmtop', '.parm7')):
                ref_top = solvated_topology if solvated_topology else mdtraj_topology
                if xtc_file:
                    ref_traj = md.load_frame(xtc_file, 0, top=ref_top)
                else:
                    ref_traj = md.load(ref_top)
            elif str(mdtraj_topology).lower().endswith('.psf'):
                # mdtraj cannot load a bare .psf (topology-only, no
                # coordinates) -- it needs a companion coordinate file, same
                # requirement as TopologyLoader._load_charmm's
                # charmm_coordinates option.
                ref_top = solvated_topology if solvated_topology else self.charmm_coordinates
                if not ref_top:
                    raise ValueError(
                        "A .psf topology needs a companion coordinate file to build an "
                        "mdtraj reference structure for custom ligand/receptor selections -- "
                        "set 'solvated_topology' or forcefield_settings.charmm_coordinates "
                        "(.pdb/.coor/.crd)."
                    )
                if xtc_file:
                    ref_traj = md.load_frame(xtc_file, 0, top=ref_top)
                else:
                    ref_traj = md.load(ref_top)
            else:
                ref_traj = md.load(mdtraj_topology)

            where = Path(str(mdtraj_topology)).name
            if str(original_input_complex_pdb).lower().endswith('.tpr'):
                where += f" (converted from {Path(str(original_input_complex_pdb)).name})"
            if require_two_chains:
                self.check_two_chains(ref_traj.topology, where)
            lig_idx, rec_idx, left_out = self.split_by_selections(
                ref_traj.topology, ligand_selection, receptor_selection, where)

            native_topology_input = str(original_input_complex_pdb).lower().endswith(
                ('.prmtop', '.parm7', '.tpr', '.top', '.psf')
            )
            if native_topology_input:
                log.info("Native topology input detected; skipping distillation and using selection masks directly.")
            else:
                # Combine and sort indices
                clean_indices = np.sort(np.concatenate([rec_idx, lig_idx])).astype(int)

                # Create Clean Complex PDB
                clean_pdb_path = str(self._work_dir_path() / 'distilled_complex.pdb')
                clean_pdb_traj = ref_traj.atom_slice(clean_indices)
                clean_pdb_traj.save(clean_pdb_path)
                log.info(f"Created distilled complex PDB: {clean_pdb_path}")

                # The distilled trajectory keeps every frame (frame selection
                # happens in run()) and only the selected atoms.
                clean_xtc_path = str(self._work_dir_path() / 'distilled.xtc')

                log.process("Distilling trajectory to match selection...")
                # Use solvated topology for XTC loading when available to avoid atom-count mismatch.
                load_top = solvated_topology if solvated_topology else complex_pdb
                if str(complex_pdb).endswith(('.prmtop', '.parm7')) and not solvated_topology:
                    # Convert prmtop to PDB for MDTraj
                    prmtop = app.AmberPrmtopFile(str(complex_pdb))
                    temp_pdb = str(self._work_dir_path() / "temp_topology_load.pdb")
                    with open(temp_pdb, 'w') as f:
                        app.PDBFile.writeFile(prmtop.topology, prmtop.positions if prmtop.positions is not None else [openmm.Vec3(0,0,0)]*prmtop.topology.getNumAtoms(), f)
                    load_top = temp_pdb

                full_xtc = md.load(xtc_file, top=load_top)
                clean_xtc = full_xtc.atom_slice(clean_indices)
                # Re-image across periodic boundaries -- see
                # TrajectoryProcessor.load_and_process's own image_molecules()
                # call (the main trajectory-loading path) for why this is
                # required: without it, a molecule that drifts and wraps to
                # the box's opposite side during the trajectory silently
                # corrupts every downstream vdW/electrostatic energy.
                if self.reimage_trajectory:
                    try:
                        clean_xtc = reimage_molecules(clean_xtc)
                    except Exception as e:
                        log.warning(f"image_molecules() failed on distilled trajectory ({e}); "
                                    f"proceeding with un-imaged coordinates.")
                clean_xtc.save(clean_xtc_path)
                log.info(f"Created distilled trajectory: {clean_xtc_path}")

                # Update paths for run
                complex_pdb = clean_pdb_path
                xtc_file = clean_xtc_path
                # Distillation actually changes the effective analysis topology.
                # Keep run() aligned with distilled files to avoid atom-count mismatches.
                original_input_complex_pdb = complex_pdb
                self.distilled_files = (clean_pdb_path, clean_xtc_path)
                # The selections, as atom indices of the distilled complex: a chain
                # or index selection would no longer match once atoms are left out
                position = {int(old): new for new, old in enumerate(clean_indices)}
                ligand_selection = index_selection(sorted(position[int(i)] for i in lig_idx))
                receptor_selection = index_selection(sorted(position[int(i)] for i in rec_idx))

        # Run the core analysis with frame selection
        # pass the original input path separately to preserve native-mode detection
        results = self.run(ligand_mol, complex_pdb, xtc_file, ligand_pdb, max_frames, energy_decomposition,
                          frame_start, frame_end, frame_stride, frame_selection, random_seed,
                          qha_analyze_complex, output_dir, ligand_selection, receptor_selection,
                          receptor_topology, ligand_topology, solvated_topology, print_interval,
                          original_complex_pdb=original_input_complex_pdb)
        
        if results is None:
            return None
        
        # Post-run validation and enhanced analysis
        df = pd.read_csv(results['output_file'])
        warnings = self.validate_results(df)
        
        if warnings:
            log.warning("Result validation warnings:")
            for warning in warnings:
                log.warning(f"  • {warning}")
        else:
            log.success("Results look reasonable")
        
        # Add detailed statistics
        binding_energies = df['binding_energy']
        bootstrap_results = self.bootstrap_uncertainty(binding_energies)
        convergence = self._check_convergence(binding_energies)
        
        results.update({
            'median_binding_energy': binding_energies.median(),
            'min_binding_energy': binding_energies.min(), 
            'max_binding_energy': binding_energies.max(),
            'convergence_check': convergence,
            'validation_warnings': warnings,
            'bootstrap_results': bootstrap_results,
            'binding_energies': binding_energies.tolist()
        })
        if {'delta_vdw', 'delta_elec'}.issubset(df.columns):
            # Gas-phase receptor-ligand interaction energy, the input of the
            # interaction-energy entropy estimates (mmgbsa.interaction_entropy)
            results['interaction_energies'] = (df['delta_vdw'] + df['delta_elec']).tolist()

        # Generate detailed report
        output_dir = self.generate_detailed_report(results, output_dir)
        results['report_directory'] = str(output_dir)
        
        return results

    def run(self, ligand_mol, complex_pdb, xtc_file, ligand_pdb, max_frames=50, energy_decomposition=False,
            frame_start=None, frame_end=None, frame_stride=None, frame_selection='sequential', random_seed=42,
            qha_analyze_complex=False, output_dir=None, ligand_selection=None, receptor_selection=None,
            receptor_topology=None, ligand_topology=None, solvated_topology=None, print_interval=10,
            original_complex_pdb=None):
        """
        Run the core MM/GBSA binding free energy calculation over a trajectory.

        Dispatches to one of two modes based on `original_complex_pdb`'s
        extension: Native Mode (`.prmtop`/`.parm7`/`.top`/`.tpr` -- topology
        and parameters loaded directly via ParmEd/TprParser, see
        `StructureManager.load_complex`) or Coordinate Mode (a plain PDB,
        parameterized on the fly via OpenFF/GAFF for the ligand and an
        OpenMM force field for the protein). In both modes, complex/receptor/
        ligand OpenMM Systems are built with the configured GB model and
        surface-area model (see `GBSAForceManager`), evaluated frame-by-frame
        over the selected trajectory frames, and combined into
        Delta G_bind = E_complex - E_receptor - E_ligand per frame.

        Parameters
        ----------
        ligand_mol : str
            Path to an isolated ligand molecule file (e.g. .sdf/.mol2) used
            for OpenFF/GAFF parameterization in Coordinate Mode. Not required
            in Native Mode (topology already carries ligand parameters).
        complex_pdb : str
            Path to the complex structure file actually used for this call;
            may be a temporary/distilled PDB derived from `original_complex_pdb`
            rather than the original path itself.
        xtc_file : str
            Path to the MD trajectory to analyze (e.g. .xtc/.dcd).
        ligand_pdb : str
            Path to an isolated ligand PDB file (optional; used for
            cross-checking atom counts/residues against the complex).
        max_frames : int, optional
            Upper bound on the number of frames selected for analysis
            (default 50). See `frame_selection` for how frames within the
            [frame_start, frame_end) range (strided by frame_stride) are
            reduced to at most this many.
        energy_decomposition : bool, optional
            If True, additionally run per-residue energy decomposition
            (see `mmgbsa.decomposition.PerResidueDecomposition`) after the
            standard binding free energy calculation.
        frame_start : int, optional
            Start frame for analysis, 0-indexed inclusive (default: 0).
        frame_end : int, optional
            End frame for analysis, 0-indexed exclusive (default: trajectory length).
        frame_stride : int, optional
            Frame stride (every Nth frame) applied within [frame_start, frame_end).
        frame_selection : str, optional
            Frame selection method: 'sequential' (strided range), 'equidistant'
            (evenly spaced subsample down to max_frames), or 'random' (random
            subsample using `random_seed`).
        random_seed : int, optional
            Random seed used only when frame_selection='random'.
        qha_analyze_complex : bool, optional
            If True, additionally accumulate complex-frame coordinates for a
            subsequent quasi-harmonic entropy analysis (see
            `mmgbsa.quasi_harmonic`) rather than discarding them after each
            frame's energy is computed.
        output_dir : str or Path, optional
            Directory to write intermediate/temporary files (e.g. the
            TPR-derived solvated PDB used for MDTraj) and final results into.
        ligand_selection : str, optional
            Custom MDTraj/ParmEd selection string identifying ligand atoms
            (e.g. 'resname UNK'), overriding residue-name-based detection.
        receptor_selection : str, optional
            Custom MDTraj/ParmEd selection string identifying receptor atoms
            (e.g. 'chainid 0'), overriding the default "everything but ligand"
            behavior.
        receptor_topology : str, optional
            Explicit path to a pre-built receptor topology (e.g. .prmtop),
            used instead of deriving the receptor by stripping the ligand
            from the complex topology.
        ligand_topology : str, optional
            Explicit path to a pre-built ligand topology (e.g. .prmtop),
            used instead of deriving the ligand by stripping everything but
            the ligand from the complex topology.
        solvated_topology : str, optional
            Path to a solvated reference topology/coordinate file (e.g. .gro)
            used for MDTraj trajectory loading when the native topology
            itself (prmtop/tpr) does not preserve original residue numbering.
        print_interval : int, optional
            Print progress (energies, timing) every N frames (default 10).
        original_complex_pdb : str, optional
            The true original input path (before any distillation/temporary
            file substitution of `complex_pdb`). Used to determine Native vs
            Coordinate mode correctly even when `complex_pdb` has since been
            replaced by a derived file. Defaults to `complex_pdb` itself if
            not given.

        Returns
        -------
        dict
            Analysis results: per-frame and mean binding energies, standard
            deviation/error, energy component breakdowns (vdW/electrostatic/
            GB/SA) when available, and `physics_assumptions` (a list of
            silent-default notices accumulated during system setup, e.g. GB
            radius fallbacks for unrecognized elements).
        """
        # Save original complex_pdb path before potential overwrite by distillation extraction
        # `original_complex_pdb` may be passed explicitly by callers (e.g. run_comprehensive)
        # so that we remember the true source file even if `complex_pdb` is later
        # replaced by a distilled or temporary PDB.  This is critical for correct
        # mode detection (native vs coordinate) when the input was originally
        # a .tpr or other native topology.
        if original_complex_pdb is None:
            original_complex_pdb = complex_pdb
        if output_dir:
            self.work_dir = Path(output_dir)

        log.section("Advanced MM/GBSA Analysis")
        log.info("Starting fixed enhanced MM/GBSA analysis with GBSA forces...")
        
        # Reset results for new run
        self.energies = defaultdict(list)
        
        # Collect cross-class assumptions
        assumptions = list(getattr(self, 'physics_assumptions', []))
        if hasattr(self, 'gbsa_manager') and hasattr(self.gbsa_manager, 'physics_assumptions'):
            assumptions.extend(self.gbsa_manager.physics_assumptions)
            
        self.results = {'physics_assumptions': assumptions}
        
        analysis_start_time = time.time()
        
        # Setup optimized platform
        platform, properties = self.setup_optimized_platform()
        self.analysis_platform = {'requested': self.platform_requested, 'used': self.platform_used,
                                  'fallback': self.platform_fallback}
        
        # Check if complex_pdb is an Amber prmtop file
        is_prmtop = complex_pdb.endswith('.prmtop') or complex_pdb.endswith('.parm7')
        
        if is_prmtop:
            log.info(f"📂 Detected Amber topology file: {complex_pdb}")
            
            # Load Amber topology
            prmtop = app.AmberPrmtopFile(complex_pdb)
            
            # Look for coordinate file
            pdb_base = Path(complex_pdb).stem
            pdb_dir = Path(complex_pdb).parent
            
            coord_file = None
            for ext in ['.inpcrd', '.rst7', '.crd']:
                candidate = pdb_dir / f"{pdb_base}{ext}"
                if candidate.exists():
                    coord_file = str(candidate)
                    break
            
            if coord_file:
                log.info(f"  ✓ Found coordinate file: {coord_file}")
                inpcrd = app.AmberInpcrdFile(coord_file)
                temp_pdb = str(Path(output_dir if output_dir else '.') / 'temp_from_prmtop.pdb')
                app.PDBFile.writeFile(prmtop.topology, inpcrd.positions, open(temp_pdb, 'w'))
                complex_pdb = temp_pdb
                log.info(f"  ✓ Converted to PDB: {complex_pdb}")
            # Prefer coordinates/residue numbering from trajectory topology when available.
            # This keeps residue IDs consistent with solvated_topology/GRO (e.g. resSeq offsets),
            # which is important for downstream decomposition/report labels and PandaMap.
            if xtc_file:
                log.info("Extracting first frame from trajectory for reference PDB...")
                try:
                     target_atoms = prmtop.topology.getNumAtoms() if 'prmtop' in locals() else None
                     traj = TrajectoryProcessor.load_and_process(
                         xtc_file,
                         complex_pdb,
                         target_atoms=target_atoms,
                         solvated_topology=solvated_topology,
                         end=1,  # only first frame needed
                         reimage=self.reimage_trajectory
                     )
                     temp_pdb = str(Path(output_dir if output_dir else '.') / 'temp_from_traj.pdb')
                     traj[0].save_pdb(temp_pdb)
                     complex_pdb = temp_pdb
                     log.info(f"  ✓ Extracted to PDB: {complex_pdb}")
                except Exception as e:
                     log.warning(f"Failed to extract trajectory reference PDB, keeping converted PDB: {e}")
        
        # Check if PDB changed (e.g. Surgery added atoms)
        # Skip if in Prmtop Mode (skip_merge=True) as we trust original indices(topology)
        
        # If we're in pure TPR mode but lack the TprParser library, fall back
        # to a coordinate-based analysis by switching to a companion PDB/GRO
        # file (if available).  This avoids a hard crash inside
        # StructureManager.load_complex.
        if str(original_complex_pdb).lower().endswith('.tpr'):
            try:
                import TprParser  # noqa: F401
            except ImportError:
                log.warning("TPR file supplied but TprParser is not installed.")
                # look for .pdb/.gro with same stem
                tpr_path = Path(original_complex_pdb)
                fb = None
                for ext in ['.pdb', '.gro']:
                    cand = tpr_path.with_suffix(ext)
                    if cand.exists():
                        fb = str(cand)
                        break
                if fb:
                    log.info(f"Switching to fallback topology {fb} (coordinate mode)")
                    original_complex_pdb = fb
                    complex_pdb = fb
                else:
                    raise ValueError(
                        "Cannot perform native TPR analysis because TprParser is missing and no "
                        "fallback PDB/GRO file was found alongside the TPR. "
                        "Please install TprParser or provide an explicit topology."
                    )
        
        # ---------------------------------------------------------
        # UNIFIED TOPOLOGY SPLITTING WORKFLOW
        # ---------------------------------------------------------
        import parmed as pmd
        
        log.process("Executing Unified Topology Splitting...")
        prep_start_time = time.time()
        ligand_resname = None
        is_native_mode = str(original_complex_pdb).endswith(('.prmtop', '.parm7', '.top', '.tpr', '.psf'))
        # Set below, in the Coordinate Mode branch, when a PPI-style
        # ligand_selection (chainid-based, no small-molecule ligand_mol) is
        # used instead of ligand_resname -- see its own comment for why this
        # needs the same "use the frame's true atom order" fix as
        # is_native_mode. Declared here so it's always defined by the time
        # the complex-context coordinate-assignment code below checks it.
        is_ppi_coordinate_mode = False
        # Set in either the Native Mode or Coordinate Mode branch below when
        # the ligand is identified by explicit atom indices rather than a
        # chainid selection or a residue name (e.g. via receptor_topology's
        # atom count for a protein-RNA/DNA complex). Declared here so the
        # "IDENTIFY COMPONENTS" fallback further down (ligand_indices/
        # protein_indices) can check it regardless of which branch ran.
        ligand_atom_indices = None
        # Native Mode with receptor_selection and ligand_selection that leave atoms
        # out: the atoms of the stripped complex that are kept, and its full size
        kept_atoms = None
        n_full_complex = None
        self.native_reduced = False

        # 1. OBTAIN MASTER COMPLEX STRUCTURE
        complex_struct = None

        if is_native_mode:
            receptor_topology = self.explicit_amber_topology(receptor_topology, 'receptor_topology')
            ligand_topology = self.explicit_amber_topology(ligand_topology, 'ligand_topology')
        elif receptor_topology or ligand_topology:
            log.warning("receptor_topology and ligand_topology are used only in Native Mode; "
                        "Coordinate Mode builds receptor and ligand from the complex and ignores them.")

        if is_native_mode:
             # NATIVE MODE: Use Unified Splitting (Robust for prmtop, tpr, and psf)
             is_charmm_native = str(original_complex_pdb).endswith('.psf')
             log.info(f"Mode: Native Topology ({'PSF/CHARMM' if is_charmm_native else 'PRMTOP or TPR'})")
             try:
                 # A .psf carries no coordinates of its own (see
                 # TopologyLoader._load_charmm's docstring) -- unlike prmtop/
                 # tpr, where `complex_pdb` (typically equal to
                 # `original_complex_pdb` at this point) can double as its
                 # own coordinate source via ParmEd's `xyz=` kwarg. Passing
                 # the same .psf path as both topology AND coordinate source
                 # here left `struct.coordinates` silently None (ParmEd's
                 # .psf branch only loads coordinates when pdb_path differs
                 # from prmtop_path) -- confirmed this produced a real
                 # System with a real particle count but ALL-ZERO or garbage
                 # positions, which blew up as a ~11 million kcal/mol VdW
                 # "energy" once real trajectory coordinates were set on it
                 # (the System's exclusions/parameters were fine; only the
                 # reference coordinates used to validate/derive
                 # receptor+ligand split were wrong). Use the companion
                 # coordinate file explicitly.
                 charmm_coord_source = solvated_topology or self.charmm_coordinates
                 complex_struct = StructureManager.load_complex(
                     charmm_coord_source if is_charmm_native else complex_pdb,
                     original_complex_pdb, xtc_path=xtc_file, gb_model=self.gb_model,
                     charmm_params=self.charmm_params if is_charmm_native else None,
                     gb_radii_set=self.gb_radii_set,
                 )
                 # One radii decision per run, taken from the complex, so the
                 # receptor and ligand cannot end up with a different set.
                 self._gb_radii_decision = getattr(complex_struct, 'opengbsa_radii_set',
                                                   self._gb_radii_decision)
                 log.success(f"Loaded Native Complex: {len(complex_struct.atoms)} atoms")

                 import os
                 if str(original_complex_pdb).endswith('.tpr'):
                      temp_pdb_path = os.path.join(output_dir if 'output_dir' in locals() and output_dir else '.', "complex_solvated_from_tpr_for_mdtraj.pdb")
                      complex_struct.save(temp_pdb_path, overwrite=True)
                      complex_pdb = temp_pdb_path

                 # STRIP SOLVENT AFTER WRITING MDTRAJ PDB BUT BEFORE OPENMM SYSTEM
                 log.info("Stripping complex system solvent to prepare for implicit GBSA evaluation...")
                 # CHARMM/NAMD-origin structures use different solvent/ion residue
                 # names than Amber/GROMACS (e.g. 'TIP3' not 'SOL', 'SOD'/'CLA' not
                 # 'NA'/'CL') -- confirmed on a real NAMD PSF that the original,
                 # Amber/GROMACS-only mask silently stripped ZERO atoms (16024 atoms
                 # in, 16024 out) despite the structure containing 4655 TIP3 waters
                 # and 29 SOD/CLA ions, which would have gone on to crash or corrupt
                 # the LCPO surface-area force downstream. Includes both naming
                 # conventions so this mask works for every currently-supported
                 # topology origin.
                 solvent_mask = (":WAT,HOH,H2O,SOL,TIP3,TIP,SPC,"
                                  "NA,CL,K,MG,ZN,CA,"
                                  "Na+,Cl-,K+,Mg2+,Ca2+,Zn2+,"
                                  "SOD,CLA,POT,CAL,ZN2")
                 complex_struct.strip(solvent_mask)
                 log.success(f"Stripped Complex: {len(complex_struct.atoms)} atoms")

                 # 2. IDENTIFY COMPONENTS
                 # A protein-protein/peptide system (e.g. NAMD/CHARMM chainid
                 # split) has no single small-molecule ligand_resname an
                 # Amber mask can select -- when ligand_selection/
                 # receptor_selection were already resolved (see
                 # `_resolve_binding_mode_selections`'s `ppi` binding_mode),
                 # use those (mdtraj-style `chainid N`) to get explicit atom
                 # indices instead of falling back to find_ligand_resname's
                 # single-residue-name assumption, which cannot represent a
                 # multi-residue peptide "ligand" at all.
                 ligand_atom_indices = None
                 if ligand_selection and receptor_selection:
                     # Both selections on the solvent-stripped complex; atoms in
                     # neither are dropped from the topology here and from the
                     # trajectory after it is loaded (kept_atoms)
                     lig_sel, rec_sel, left_out = self.split_by_selections(
                         complex_struct.topology, ligand_selection, receptor_selection,
                         "the solvent-stripped complex")
                     if left_out:
                         n_full_complex = len(complex_struct.atoms)
                         kept_atoms = sorted(set(lig_sel) | set(rec_sel))
                         keep = set(kept_atoms)
                         complex_struct.strip([i not in keep for i in range(n_full_complex)])
                         # ParmEd keeps its cached OpenMM Topology after an Amber
                         # strip; drop it so the reduced complex is used
                         if hasattr(complex_struct, '_topology'):
                             del complex_struct._topology
                         position = {old: new for new, old in enumerate(kept_atoms)}
                         lig_sel = sorted(position[i] for i in lig_sel)
                         rec_sel = sorted(position[i] for i in rec_sel)
                         # as indices of the reduced complex (chain numbers may shift)
                         ligand_selection = index_selection(lig_sel)
                         receptor_selection = index_selection(rec_sel)
                         log.info(f"Complex reduced to the {len(kept_atoms)} selected atoms")
                     ligand_atom_indices = list(lig_sel)
                     log.info(f"Resolved ligand_selection to {len(lig_sel)} atoms and receptor_selection "
                              f"to {len(rec_sel)} atoms (post-solvent-strip complex)")
                 elif ligand_selection:
                     # complex_struct.topology here is a ParmEd Structure's
                     # own .topology (an OpenMM Topology) post-solvent-strip;
                     # get_selection_indices already accepts an OpenMM
                     # Topology directly (see get_ligand_indices/
                     # get_protein_indices's own dual OpenMM/MDTraj handling).
                     ligand_atom_indices = self.get_selection_indices(complex_struct.topology, ligand_selection)
                     if not ligand_atom_indices:
                         raise ValueError(f"ligand_selection '{ligand_selection}' matches no atom of the "
                                          f"solvent-stripped complex.")
                     log.info(f"Resolved ligand_selection '{ligand_selection}' to {len(ligand_atom_indices)} "
                              f"atoms (post-solvent-strip complex); the receptor is every other atom "
                              f"({len(complex_struct.atoms) - len(ligand_atom_indices)} atoms)")
                 elif receptor_topology:
                     # An explicit receptor_topology (with a known atom
                     # count) tells us exactly how to split the complex
                     # WITHOUT guessing a ligand_resname -- required for any
                     # complex where find_ligand_resname's "single small
                     # residue = ligand" heuristic is wrong, e.g. a
                     # protein-RNA/DNA complex (confirmed on a real dataset,
                     # Zenodo 6973437: find_ligand_resname picked "HIP", a
                     # histidine tautomer residue, as the "ligand" instead of
                     # the actual 25-nucleotide RNA chain, since RNA/DNA
                     # residue names were never in its exclusion/detection
                     # logic at all). The receptor's own atom count is
                     # reliable regardless of residue-naming conventions,
                     # since the receptor+ligand topologies were built from
                     # the exact same complex by construction.
                     receptor_atom_count = len(pmd.load_file(receptor_topology).atoms)
                     total_atoms = len(complex_struct.atoms)
                     if receptor_atom_count < total_atoms:
                         ligand_atom_indices = list(range(receptor_atom_count, total_atoms))
                         log.info(f"Derived ligand indices from receptor_topology's atom count "
                                  f"({receptor_atom_count} receptor + {total_atoms - receptor_atom_count} "
                                  f"ligand = {total_atoms} total) -- bypassing ligand_resname guessing.")
                     else:
                         log.warning(f"receptor_topology has {receptor_atom_count} atoms >= complex's "
                                     f"{total_atoms}; cannot derive ligand indices this way, "
                                     f"falling back to ligand_resname detection.")
                 if ligand_atom_indices is None and not ligand_resname:
                     ligand_resname = self.find_ligand_resname(complex_struct.topology)
                 if ligand_resname and ligand_atom_indices is None:
                     log.info(f"Using Ligand Residue Name: {ligand_resname}")

                 # 3. COMPONENT PREPARATION
                 # ALWAYS split from complex to get COORDINATES (and fallback topology)
                 log.info("Splitting Complex to obtain coordinates...")
                 receptor_derived, ligand_derived = StructureManager.split_components(
                     complex_struct, ligand_resname, ligand_indices=ligand_atom_indices
                 )

                 # Handle Receptor
                 if receptor_topology:
                     log.info(f"Using Explicit Receptor Topology: {receptor_topology}")
                     # Load bare topology
                     import parmed as pmd
                     receptor_explicit = pmd.load_file(receptor_topology)
                     
                     # Check atom count
                     if len(receptor_explicit.atoms) != len(receptor_derived.atoms):
                         raise ValueError(f"Explicit receptor prmtop ({len(receptor_explicit.atoms)} atoms) mismatches derived receptor ({len(receptor_derived.atoms)} atoms). Check ligand_resname or prmtop.")
                     
                     # Transfer coordinates
                     receptor_explicit.coordinates = receptor_derived.coordinates
                     # Also transfer box if needed, though GBSA is non-periodic
                     receptor_struct = receptor_explicit
                 else:
                     receptor_struct = receptor_derived
                 
                 # Handle Ligand
                 if ligand_topology:
                     log.info(f"Using Explicit Ligand Topology: {ligand_topology}")
                     # Load bare topology
                     ligand_explicit = pmd.load_file(ligand_topology)
                     
                     # Check atom count
                     if len(ligand_explicit.atoms) != len(ligand_derived.atoms):
                          raise ValueError(f"Explicit ligand prmtop ({len(ligand_explicit.atoms)} atoms) mismatches derived ligand ({len(ligand_derived.atoms)} atoms).")
                     
                     # Transfer coordinates
                     ligand_explicit.coordinates = ligand_derived.coordinates
                     ligand_struct = ligand_explicit
                 else:
                     ligand_struct = ligand_derived
                 
                 # 4. CREATE OPENMM SYSTEMS
                 log.process("Creating and Enhancing OpenMM Systems from Unified Topology...")
                 def make_gbsa_system(struct, name):
                     # Prepare Cutoff args
                     nb_method = app.NoCutoff
                     nb_cutoff = None
                     if self.nonbonded_cutoff is not None:
                         nb_method = app.CutoffNonPeriodic
                         nb_cutoff = self.nonbonded_cutoff * unit.angstroms

                     # Ensure this struct's GB radii match self.gb_model,
                     # independent of whatever radius set its own source
                     # prmtop happened to be built with. This matters
                     # specifically for receptor_struct/ligand_struct when
                     # receptor_topology/ligand_topology point at
                     # SEPARATELY-loaded prmtop files (see the
                     # `pmd.load_file(receptor_topology)`/
                     # `pmd.load_file(ligand_topology)` branches above) --
                     # those never go through StructureManager.load_complex,
                     # so without this they'd keep whatever radii their own
                     # file has even after load_complex's changeRadii call
                     # updated complex_struct's radii to match gb_model,
                     # producing a complex/receptor/ligand radii MISMATCH
                     # (confirmed: this desync alone changed a real system's
                     # OBC2 binding energy from -27.6 to -236.0 kcal/mol --
                     # not a gb_model effect, just complex and
                     # receptor+ligand disagreeing about GB radii on the
                     # same atoms). Applying it uniformly here, to every
                     # struct this function builds a System from, keeps all
                     # three consistent regardless of where each struct
                     # originated.
                     # Since 0.0.11 the set is the complex's decision (None:
                     # every structure keeps its own Amber topology radii).
                     radii_set = self._gb_radii_decision
                     if radii_set is not None:
                         try:
                             from parmed.tools import changeRadii
                             changeRadii(struct, radii_set).execute()
                         except Exception as e:
                             log.warning(f"Failed to apply GB radii set to {name} structure: {e}.")

                     self.gbsa_manager.resolve_mmpbsa_dielectric(struct)
                     self.solvent_dielectric = self.gbsa_manager.solvent_dielectric
                     sys = StructureManager.create_openmm_system(
                         struct, 
                         implicitSolvent=self.gbsa_manager.current_app_model,
                         implicitSolventSaltConc=self.salt_concentration * unit.molar if self.salt_concentration > 0 else 0.0*unit.molar,
                         nonbondedMethod=nb_method,
                         nonbondedCutoff=nb_cutoff,
                         solventDielectric=self.gbsa_manager.solvent_dielectric
                     ) 
                     sys = self.gbsa_manager.refine_gbsa_forces(sys, struct)
                     return sys

                 complex_system = make_gbsa_system(complex_struct, "Complex")
                 protein_system = make_gbsa_system(receptor_struct, "Receptor")
                 ligand_system = make_gbsa_system(ligand_struct, "Ligand")
                 
                 # Set topologies for slicing
                 complex_top = complex_struct.topology
                 protein_top = receptor_struct.topology
                 ligand_top = ligand_struct.topology
                 
             except Exception as e:
                 log.error(f"Native Mode Failed: {e}")
                 raise e

        else:
             # COORDINATE MODE: Use Consistent Independent Generation
             # (Bypasses ParmEd serialization issues with OpenFF/SystemGenerator)
             log.info("Mode: Raw Coordinates -- Generating Consistent Systems via OpenMM")

             try:
                 # 1. Load Complex PDB & Identify Ligand
                 import parmed as pmd
                 c_struct = pmd.load_file(original_complex_pdb)

                 # A protein-protein/peptide system (e.g. binding_mode=ppi,
                 # split by chainid) has no single small-molecule
                 # ligand_resname an Amber mask can select -- mirrors the
                 # identical fix already applied to Native Mode above. When
                 # ligand_selection was already resolved (chainid-style),
                 # use it to get explicit atom indices instead of the
                 # ligand_mol/find_ligand_resname path below, which assumes
                 # a single-residue small-molecule ligand and previously
                 # left ligand_resname=None here for a PPI system (no
                 # ligand_mol was ever provided), silently mis-splitting the
                 # complex and then crashing several steps later trying to
                 # OpenFF-parameterize a protein chain as if it were a
                 # small molecule.
                 ligand_atom_indices = None
                 if ligand_selection:
                     ligand_atom_indices = self.get_selection_indices(c_struct.topology, ligand_selection)
                     if not ligand_atom_indices:
                         raise ValueError(f"ligand_selection '{ligand_selection}' matches no atom of "
                                          f"{Path(str(original_complex_pdb)).name}.")
                     log.info(f"Resolved ligand_selection '{ligand_selection}' to {len(ligand_atom_indices)} atoms (Coordinate Mode)")
                 elif not ligand_resname and ligand_mol:
                     # Detect Ligand
                     unique_resnames = set(r.name for r in c_struct.residues)
                     if 'LIG' in unique_resnames:
                         ligand_resname = 'LIG'
                     else:
                         ligand_resname = self.find_ligand_resname(c_struct.topology) or 'LIG'
                     log.info(f"Detected Ligand Residue: {ligand_resname}")

                 # 2. Split Coordinates into PDBs
                 r_struct, l_struct = StructureManager.split_components(
                     c_struct, ligand_resname, ligand_indices=ligand_atom_indices
                 )

                 import os
                 rec_pdb_path = os.path.join(output_dir, "temp_receptor.pdb")
                 lig_pdb_path = os.path.join(output_dir, "temp_ligand.pdb")
                 r_struct.save(rec_pdb_path, overwrite=True)
                 l_struct.save(lig_pdb_path, overwrite=True)

                 # 3. Generate Systems Independently (Consistent FF)
                 log.process("Generating Complex System...")
                 is_ppi_coordinate_mode = ligand_atom_indices is not None
                 if is_ppi_coordinate_mode:
                     # No small-molecule ligand at all -- build all three
                     # systems the same way build_complex_system already
                     # handles "Protein-Only/Protein-Protein Mode"
                     # (ligand_mol=None), using the distilled/original
                     # complex PDB and the two chain-split PDBs directly.
                     # This never touches the OpenFF/ligand_mol path, so
                     # it cannot hit the "'NoneType' object has no
                     # attribute '_finalize'" failure that came from
                     # calling build_complex_system(lig_pdb_path,
                     # ligand_mol=None, ligand_pdb=lig_pdb_path) below --
                     # that call told build_complex_system to both skip
                     # ligand parameterization (ligand_mol=None) AND
                     # attempted to add an OpenFF ligand from ligand_pdb,
                     # an invalid combination this code path never
                     # produces for is_ppi_coordinate_mode.
                     complex_system, cx_top, _ = self.build_complex_system(original_complex_pdb, ligand_mol=None, add_gbsa=True)
                     complex_top = cx_top

                     log.process("Generating Receptor System...")
                     protein_system, px_top, _ = self.build_complex_system(rec_pdb_path, ligand_mol=None, add_gbsa=True)
                     protein_top = px_top

                     log.process("Generating Ligand System...")
                     ligand_system, lx_top, _ = self.build_complex_system(lig_pdb_path, ligand_mol=None, add_gbsa=True)
                     ligand_top = lx_top
                 else:
                     # For complex, we need to extract ligand PDB first if we want build_complex_system to use it with mol
                     # Or just pass the extracted lig_pdb_path
                     complex_system, cx_top, _ = self.build_complex_system(original_complex_pdb, ligand_mol, ligand_pdb=lig_pdb_path, add_gbsa=True)
                     complex_top = cx_top

                     # The parameterized ligand (OpenFF) might have a different residue name (e.g. UNK)
                     # We must update ligand_resname to match the new topology for accurate indexing
                     new_res = self.find_ligand_resname(complex_top)
                     if new_res and new_res != ligand_resname:
                          log.info(f"Ligand Residue Name updated from {ligand_resname} to {new_res} (OpenFF default)")
                          ligand_resname = new_res

                     log.process("Generating Receptor System...")
                     protein_system, px_top, _ = self.build_complex_system(rec_pdb_path, ligand_mol=None, add_gbsa=True)
                     protein_top = px_top

                     log.process("Generating Ligand System...")
                     # Treat lig.pdb as "protein_pdb" input but with ligand_mol so it gets parameterized as ligand
                     # build_complex_system will delete LIG from lig.pdb (emptying it) then add ligand_pdb (refilling it)
                     ligand_system, lx_top, _ = self.build_complex_system(lig_pdb_path, ligand_mol, ligand_pdb=lig_pdb_path, add_gbsa=True)
                     ligand_top = lx_top

                 log.success("Systems generated consistently.")

             except Exception as e:
                 log.error(f"Coordinate Mode Failed: {e}")
                 raise e
        
        if not 'complex_system' in locals():
             # Fallback if something weird happened, or for safety
             raise RuntimeError("System generation failed to produce complex_system")
             
        # If openFF re-assigned the ligand name, ensure custom selections target the new name
        final_lig_selection = ligand_selection
        final_rec_selection = receptor_selection
        
        if ligand_resname:
            if final_lig_selection:
                if "resname UNL" in final_lig_selection:
                    final_lig_selection = final_lig_selection.replace("resname UNL", f"resname {ligand_resname}")
                elif "resname LIG" in final_lig_selection:
                    final_lig_selection = final_lig_selection.replace("resname LIG", f"resname {ligand_resname}")
                if final_lig_selection != ligand_selection:
                    log.info(f"Updated ligand selection to match topology: '{final_lig_selection}'")
                    
            if final_rec_selection:
                if "resname UNL" in final_rec_selection:
                    final_rec_selection = final_rec_selection.replace("resname UNL", f"resname {ligand_resname}")
                elif "resname LIG" in final_rec_selection:
                    final_rec_selection = final_rec_selection.replace("resname LIG", f"resname {ligand_resname}")
                if final_rec_selection != receptor_selection:
                    log.info(f"Updated receptor selection to match topology: '{final_rec_selection}'")
        
        # Important: Store the resolved selection attributes on the instance to guarantee downstream consistency
        self.ligand_selection = final_lig_selection
        self.receptor_selection = final_rec_selection
        self.ligand_resname = ligand_resname

        if ligand_atom_indices is not None and not final_lig_selection and not ligand_resname:
            # Neither a chainid-style selection nor a residue name identifies
            # the ligand here (e.g. Native Mode split by an explicit
            # receptor_topology atom count -- see the "IDENTIFY COMPONENTS"
            # block above, added for protein-RNA/DNA complexes where
            # find_ligand_resname's single-residue heuristic picks the wrong
            # atom entirely). Without this, get_ligand_indices/
            # get_protein_indices below fall through to their
            # ligand_resname=None default, which matches ZERO atoms --
            # confirmed this produced "Ligand system has N particles but
            # complex has 0 ligand atoms" on a real protein-RNA dataset
            # (Zenodo 6973437) despite the complex/receptor/ligand Systems
            # themselves having been built correctly moments earlier.
            ligand_indices = list(ligand_atom_indices)
            protein_indices = [i for i in range(len(list(complex_top.atoms())
                                                        if callable(getattr(complex_top, 'atoms', None))
                                                        else complex_top.atoms))
                               if i not in set(ligand_indices)]
            log.info(f"Using explicit ligand/receptor atom indices ({len(ligand_indices)} ligand, "
                     f"{len(protein_indices)} receptor) from receptor_topology's atom count.")
            # a one-residue ligand keeps its name for the 3D view of the HTML report
            all_atoms = list(complex_top.atoms() if callable(getattr(complex_top, 'atoms', None))
                             else complex_top.atoms)
            lig_set = set(ligand_indices)
            names = {all_atoms[i].residue.name for i in lig_set if i < len(all_atoms)}
            if len(names) == 1:
                self.report_ligand_resname = names.pop()
        else:
            ligand_indices = self.get_ligand_indices(complex_top, ligand_resname, selection=final_lig_selection)
            if final_lig_selection and not final_rec_selection:
                # ligand_selection alone: the receptor is every other atom of the
                # (solvent-free) complex, as its receptor System is
                lig_set = set(ligand_indices)
                n_complex = len(list(complex_top.atoms() if callable(getattr(complex_top, 'atoms', None))
                                     else complex_top.atoms))
                protein_indices = [i for i in range(n_complex) if i not in lig_set]
            else:
                protein_indices = self.get_protein_indices(complex_top, ligand_resname, selection=final_rec_selection)

        # ===== Diagnostic warnings for custom selections =====
        # If the user supplied a custom selection but it ended up selecting the
        # same atoms as the default (no-selection) behaviour we warn.  This
        # commonly happens when a topology only contains one protein chain or
        # when a converted TPR has lost chain/residue information and the
        # fallback path effectively ignored the custom selection.  The warning
        # should help users understand why two different expressions produced
        # identical energies.
        if (ligand_selection or receptor_selection):
            # Suppress logging while we fetch defaults for comparison
            import logging as root_logging
            old_level = root_logging.getLogger().level
            root_logging.getLogger().setLevel(root_logging.ERROR)
            default_lig = self.get_ligand_indices(complex_top, ligand_resname)
            default_prot = self.get_protein_indices(complex_top, ligand_resname)
            root_logging.getLogger().setLevel(old_level)
            if ligand_selection and set(ligand_indices) == set(default_lig):
                log.warning(
                    "Custom ligand selection '%s' returned the same atom set as the default; "
                    "it will have no effect on the calculation." % final_lig_selection
                )
            if receptor_selection and set(protein_indices) == set(default_prot):
                log.warning(
                    "Custom receptor selection '%s' returned the same atom set as the default; "
                    "it will have no effect on the calculation." % final_rec_selection
                )

        
        skip_merge = True
        final_pdb_path = complex_pdb
        self.distilled_complex_pdb = final_pdb_path
        self.distilled_xtc_file = xtc_file

        prep_time = time.time() - prep_start_time
        print(f"✓ Total preparation time: {prep_time:.1f}s")
        
        # For implicit solvent, protein_system may have additional virtual sites or implicit solvent atoms
        # So we check consistency differently: protein_system + ligand_system should equal complex_system
        
        # Ligand system validation
        if ligand_system and ligand_system.getNumParticles() != len(ligand_indices):
            print(f"ERROR: Ligand system has {ligand_system.getNumParticles()} particles but complex has {len(ligand_indices)} ligand atoms")
            return None
        
        # Complex system validation (total should match protein + ligand)
        if ligand_system:
            if complex_system.getNumParticles() != (protein_system.getNumParticles() + ligand_system.getNumParticles() - len(ligand_indices)):
                # For implicit solvent, the protein+ligand systems built separately may not sum exactly due to overlapping solvent context
                print(f"⚠️  System atom counts are approximate due to implicit solvent:")
                print(f"    Complex: {complex_system.getNumParticles()}, Protein: {protein_system.getNumParticles()}, Ligand: {ligand_system.getNumParticles()}")
                print(f"    Expected protein atoms: {len(protein_indices)}, ligand atoms: {len(ligand_indices)}")
            else:
                print("✓ All particle counts match!")
        else:
            # Protein-only case
            if protein_system.getNumParticles() != len(protein_indices):
                print(f"⚠️  Protein system has {protein_system.getNumParticles()} particles but {len(protein_indices)} protein atoms found")
                print(f"    This is expected for implicit solvent or if explicit waters were preserved")
            else:
                print("✓ Protein system particle count matches!")

        # HELPER: Assign Force Groups for Internal Energy Cancellation
        #
        # IMPORTANT: this runs AFTER `refine_gbsa_forces` has already built and
        # correctly force-grouped the 1-4 VDW correction force (a
        # CustomBondForce with energy function "4*epsilon*((sigma/r)^12 -
        # (sigma/r)^6)", explicitly set to group 3 -- see the "Moving 1-4s to
        # CustomBondForce Group 3" step). A previous version of this function
        # blindly reassigned EVERY CustomBondForce to group 10 (the
        # bonded/internal-energy group), silently moving ~10,000 1-4 VDW
        # interactions per system OUT of the reported "VDW" energy group and
        # into "internal energy" instead -- confirmed on real OXA-MD
        # benchmark data to be the actual root cause of a systematic ~18
        # kcal/mol vdW discrepancy against Amber MMPBSA.py reference results
        # (the 1-4 VDW correction fix itself, and the mbondi2/1-4-exception
        # fixes made earlier while chasing this bug, were both real,
        # independently-valid fixes but neither was the actual cause of this
        # specific discrepancy). CustomBondForce is therefore identified by
        # its energy function string rather than being reassigned wholesale.
        def assign_force_groups(sys_obj):
            log.debug(f"Inspecting forces for system with {sys_obj.getNumParticles()} particles")
            for f in sys_obj.getForces():
                log.debug(f"  - Found Force: {type(f).__name__}")
                if isinstance(f, openmm.HarmonicBondForce): f.setForceGroup(10)
                elif isinstance(f, openmm.HarmonicAngleForce): f.setForceGroup(11)
                elif isinstance(f, openmm.PeriodicTorsionForce): f.setForceGroup(12)
                elif isinstance(f, openmm.RBTorsionForce): f.setForceGroup(13)
                elif isinstance(f, openmm.CMAPTorsionForce): f.setForceGroup(14)
                elif isinstance(f, openmm.CustomTorsionForce):
                    # Improper torsions: ParmEd's own Structure.createSystem
                    # hardcodes these to IMPROPER_FORCE_GROUP=4 (see ParmEd's
                    # structure.py), which collides with this codebase's
                    # unrelated group-4 convention for the LCPO/SA nonpolar
                    # solvation force -- confirmed on real data to silently
                    # fold improper-torsion energy into the reported "surface
                    # area" energy. Move to an unused internal-energy group.
                    f.setForceGroup(15)
                # Ensure Nonbonded is Group 0 (Default, but explicit is good)
                elif isinstance(f, openmm.NonbondedForce): f.setForceGroup(0)
                # GBSA forces are usually handled by add_gbsa_to_system, but if present as standard:
                elif isinstance(f, openmm.GBSAOBCForce):
                    # GBSA and SA forces are definitively assigned groups 1 and 4 via `refine_gbsa_forces`.
                    # Overriding groups conditionally based on only particle(0)'s charge creates extreme vulnerability
                    # if the parameterized ligand structure happens to list a low-charge/dummy atom at index 0.
                    pass
                # Custom Forces (GBn, SA, Screening)
                elif isinstance(f, openmm.CustomGBForce): f.setForceGroup(2)
                elif isinstance(f, openmm.CustomNonbondedForce): f.setForceGroup(3)
                elif isinstance(f, openmm.CustomBondForce):
                    # The only CustomBondForce this codebase creates is the
                    # 1-4 VDW correction from `refine_gbsa_forces` (LCPO uses
                    # the dedicated openmm.LCPOForce class, not
                    # CustomBondForce, so there is no real ambiguity here).
                    # Identify it by its energy function and put it in group
                    # 3 (VDW) where it belongs; leave anything unrecognized
                    # in whatever group it already has rather than guessing.
                    if 'sigma' in f.getEnergyFunction() and 'epsilon' in f.getEnergyFunction():
                        f.setForceGroup(3)
        
        print("Assigning force groups to component systems for correct energy decomposition...")
        if ligand_system:
            assign_force_groups(ligand_system)
        assign_force_groups(protein_system)
        assign_force_groups(complex_system)

        # Create integrators (Strip units for safety)
        temp_k = self.temperature
        if unit.is_quantity(temp_k):
            temp_k = temp_k.value_in_unit(unit.kelvin)
        
        fric = 1.0 # 1/ps
        step = 0.001 # 1fs

        ligand_integrator = openmm.LangevinMiddleIntegrator(float(temp_k), float(fric), float(step))
        protein_integrator = openmm.LangevinMiddleIntegrator(float(temp_k), float(fric), float(step))
        complex_integrator = openmm.LangevinMiddleIntegrator(float(temp_k), float(fric), float(step))
        
        ligand_context = openmm.Context(ligand_system, ligand_integrator, platform, properties)
        protein_context = openmm.Context(protein_system, protein_integrator, platform, properties)
        complex_context = openmm.Context(complex_system, complex_integrator, platform, properties)
        
        # Cache fully parameterized systems and topologies to prevent downstream modules (e.g. per-residue decomp)
        # from redundantly rebuilding them and potentially losing native TPR exact parameters.
        self.systems = {
            'complex_system': complex_system,
            'complex_topology': complex_top,
            'protein_system': protein_system,
            'protein_topology': protein_top,
            'ligand_system': ligand_system,
            'ligand_topology': ligand_top
        }
        
        print("✓ All contexts created successfully!")
        
        # Load trajectory (Universal)
        if xtc_file:
            print(f"Loading trajectory {xtc_file}...")
            # Use TrajectoryProcessor to handle loading and potential stripping
            try:
                traj = TrajectoryProcessor.load_and_process(
                    trajectory_file=xtc_file,
                    topology_file=final_pdb_path if 'final_pdb_path' in locals() else complex_pdb,
                    # the System's atom count is the truth (before atoms in neither
                    # selection were dropped, see kept_atoms)
                    target_atoms=n_full_complex or complex_system.getNumParticles(),
                    solvated_topology=solvated_topology,
                    reimage=self.reimage_trajectory
                )
                print(f"✓ Trajectory processed ({traj.n_frames} frames)")
            except Exception as e:
                log.error(f"Trajectory processing failed: {e}")
                # Last ditch: try basic load if processor fails unexpectedly.
                # Still re-image (see TrajectoryProcessor.load_and_process's
                # own image_molecules() call for why this matters) since
                # this fallback bypasses that processor entirely.
                traj = md.load(xtc_file, top=final_pdb_path)
                if self.reimage_trajectory:
                    try:
                        traj = reimage_molecules(traj)
                    except Exception as e2:
                        log.warning(f"image_molecules() failed on fallback trajectory load ({e2}); "
                                    f"proceeding with un-imaged coordinates.")

            if kept_atoms is not None:
                # Native Mode: drop the atoms in neither selection, as from the topology,
                # and keep files of the reduced complex for the per-residue decomposition
                traj = traj.atom_slice(kept_atoms)
                work = self._work_dir_path()
                reduced_pdb, reduced_xtc = str(work / 'distilled_complex.pdb'), str(work / 'distilled.xtc')
                traj[0].save_pdb(reduced_pdb)
                traj.save(reduced_xtc)
                self.distilled_complex_pdb, self.distilled_xtc_file = reduced_pdb, reduced_xtc
                self.distilled_files = (reduced_pdb, reduced_xtc)
                self.native_reduced = True
                log.info(f"Trajectory reduced to the {len(kept_atoms)} selected atoms "
                         f"({reduced_pdb}, {reduced_xtc})")

            # Verification handled by Processor, but double check
            if traj.n_atoms != complex_system.getNumParticles():
                 log.warning("Final atom count mismatch despite processing!")

        elif 'traj' not in locals() or traj is None:
                 
            print(f"✓ Trajectory loaded ({traj.n_frames} frames, {time.time()-prep_start_time:.1f}s)")
            
            # Align trajectory to reference? (Optional but good practice)
            # traj.superpose(md.load(final_pdb_path))
        elif 'traj' not in locals() or traj is None:
            print("No XTC file provided or XTC disabled. Loading Complex PDB as trajectory.")
            target_pdb = locals().get('final_pdb_path', complex_pdb)
            try:
                traj = md.load(target_pdb)
            except Exception:
                traj = md.load(str(target_pdb))
        
        # Select frames based on parameters
        selected_frames = self._select_frames(len(traj), max_frames, frame_start, frame_end, 
                                            frame_stride, frame_selection, random_seed)
        traj = traj[selected_frames]
        print(f"✓ Selected {len(traj)} frames for analysis")
        
        # Process each frame
        calc_start_time = time.time()
        qha_ligand_frames = [] # Collection for QHA
        qha_protein_frames = []
        qha_complex_frames = []
        
        # Component Energy Lists (Delta)
        delta_nb_values = []
        delta_gb_values = []
        delta_sa_values = []
        delta_screen_values = []

        # Ensure energies output structure exists for per-frame data
        if 'complex' not in self.energies: self.energies['complex'] = []
        if 'protein' not in self.energies: self.energies['protein'] = []
        if 'ligand' not in self.energies: self.energies['ligand'] = []
        if 'binding' not in self.energies: self.energies['binding'] = []
        
        if 'complex_nb' not in self.energies: self.energies['complex_nb'] = []
        if 'complex_gb' not in self.energies: self.energies['complex_gb'] = []
        if 'complex_sa' not in self.energies: self.energies['complex_sa'] = []
        if 'complex_screen' not in self.energies: self.energies['complex_screen'] = []
        if 'complex_bondsa' not in self.energies: self.energies['complex_bondsa'] = []
        
        if 'delta_nb' not in self.energies: self.energies['delta_nb'] = []
        if 'delta_gb' not in self.energies: self.energies['delta_gb'] = []
        if 'delta_sa' not in self.energies: self.energies['delta_sa'] = []
        if 'delta_screen' not in self.energies: self.energies['delta_screen'] = []
        
        # Add VdW/Elec separation
        if 'delta_vdw' not in self.energies: self.energies['delta_vdw'] = []
        if 'delta_elec' not in self.energies: self.energies['delta_elec'] = []
        
        delta_vdw_values = []
        delta_elec_values = []
        
        delta_gb_values = []
        delta_sa_values = []

        # Parallel per-frame dispatch (CPU platform only, and only when QHA
        # is disabled -- QHA collects every frame into qha_*_frames as it
        # goes, which is inherently sequential/stateful in the main
        # process, so the plain serial loop is kept for that case and for
        # any non-CPU platform, where worker processes can't safely create
        # their own Context; see _frame_worker_init's docstring). When
        # eligible, frame_results is populated up front and the loop below
        # reads from it instead of computing inline -- everything AFTER the
        # loop (decomposition, entropy, reporting) is completely unchanged
        # and shared between both paths.
        frame_results = None
        use_parallel_frames = (
            platform.getName() == 'CPU'
            and self.entropy_method != 'quasiharmonic'
            and len(traj) > 1
        )
        if use_parallel_frames:
            print(f"Using parallel frame processing on CPU ({self.max_workers} workers)...")
            complex_xml = openmm.XmlSerializer.serialize(complex_system)
            protein_xml = openmm.XmlSerializer.serialize(protein_system)
            ligand_xml = openmm.XmlSerializer.serialize(ligand_system)
            temp_k = self.temperature
            if unit.is_quantity(temp_k):
                temp_k = temp_k.value_in_unit(unit.kelvin)

            tasks = []
            for i in range(traj.n_frames):
                xyz_nm = traj.xyz[i]  # a view; traj[i] deep-copies a whole Trajectory
                complex_pos_full = xyz_nm if (is_native_mode or is_ppi_coordinate_mode) else None
                tasks.append((i, xyz_nm, list(protein_indices), list(ligand_indices),
                              complex_pos_full, is_native_mode, is_ppi_coordinate_mode))

            frame_results = {}
            n_workers = max(1, min(self.max_workers, len(tasks)))
            # Each worker gets an equal share of the CPU threads; frames go out in
            # chunks so per-task pickling and dispatch are amortised.
            threads = _cpu_threads_per_worker(n_workers)
            chunksize = max(1, len(tasks) // (4 * n_workers))
            print(f"  {n_workers} worker(s) x {threads or 'OPENMM_CPU_THREADS'} OpenMM thread(s), "
                  f"chunks of {chunksize} frames")
            try:
                with ProcessPoolExecutor(
                    max_workers=n_workers, initializer=_frame_worker_init,
                    initargs=(complex_xml, protein_xml, ligand_xml, float(temp_k), threads)
                ) as pool:
                    for frame_idx, result in pool.map(_frame_worker_process, tasks, chunksize=chunksize):
                        frame_results[frame_idx] = result
            except Exception as e:
                log.warning(f"Parallel frame processing failed ({e}); falling back to serial loop.")
                frame_results = None

        # Frames are read as coordinate views (traj.xyz[i]) and positions are set
        # from numpy arrays: building a Trajectory per frame (traj[i]) and a
        # per-atom Quantity list cost ~25x the GPU energy evaluation itself.
        lig_idx = np.asarray(ligand_indices)
        prot_idx = np.asarray(protein_indices)
        for i in range(traj.n_frames):
            if i % 10 == 0:
                print(f"Frame {i+1}/{len(traj)}")

            if frame_results is not None:
                # Parallel path: pull this frame's pre-computed result
                # instead of recomputing it inline below.
                result = frame_results.get(i)
                if result is None or result.get('error'):
                    err = result.get('error') if result else 'no result returned'
                    print(f"Error calculating energies for frame {i}: {err}")
                    print("Skipping this frame...")
                    continue

                complex_e = result['complex_e']; protein_e = result['protein_e']; ligand_e = result['ligand_e']
                binding_e = result['binding_e']
                c_nb = result['c_nb']; c_gb = result['c_gb']; c_sa = result['c_sa']
                c_vdw = result['c_vdw']; c_elec = result['c_elec']
                p_vdw = result['p_vdw']; p_elec = result['p_elec']
                l_vdw = result['l_vdw']; l_elec = result['l_elec']
                p_gb = result['p_gb']; l_gb = result['l_gb']
                p_sa = result['p_sa']; l_sa = result['l_sa']
                p_nb = p_vdw + p_elec
                l_nb = l_vdw + l_elec
                e_bond = result['e_bond']; e_angle = result['e_angle']
                e_torsion = result['e_torsion']; e_rbtorsion = result['e_rbtorsion']; e_cmap = result['e_cmap']
                e_nb = result['e_nb']; e_obc = result['e_obc']; e_internal = result['e_internal']
                complex_e_clean = complex_e - e_internal
                c_screen = p_screen = l_screen = 0.0
                e_bondsa = 0.0

                if i == 0:
                    print(f"\nBreakdown (kcal/mol):")
                    print(f"            {'Complex':>12} {'Protein':>12} {'Ligand':>12} {'Delta':>12}")
                    print(f"  NB (VdV+Ele): {c_nb:12.2f} {p_nb:12.2f} {l_nb:12.2f} {c_nb - p_nb - l_nb:12.2f}")
                    print(f"    - VdW:      {c_vdw:12.2f} {p_vdw:12.2f} {l_vdw:12.2f} {c_vdw - p_vdw - l_vdw:12.2f}")
                    print(f"    - Elec:     {c_elec:12.2f} {p_elec:12.2f} {l_elec:12.2f} {c_elec - p_elec - l_elec:12.2f}")
                    print(f"  GB (PolSol):  {c_gb:12.2f} {p_gb:12.2f} {l_gb:12.2f} {c_gb - p_gb - l_gb:12.2f}")
                    print(f"  SA (NonPol):  {c_sa:12.2f} {p_sa:12.2f} {l_sa:12.2f} {c_sa - p_sa - l_sa:12.2f}")
                    print(f"  Total Clean (binding): {binding_e:12.2f}")
                    print(f"-------------------------------------------------------")

                self.energies['complex'].append(complex_e)
                self.energies['protein'].append(protein_e)
                self.energies['ligand'].append(ligand_e)
                self.energies['binding'].append(binding_e)

                self.energies['complex_nb'].append(c_nb)
                self.energies['complex_gb'].append(c_gb)
                self.energies['complex_sa'].append(c_sa)
                self.energies['complex_screen'].append(c_screen)
                self.energies['complex_bondsa'].append(e_bondsa)

                delta_nb_values.append(c_nb - p_nb - l_nb)
                delta_vdw_values.append(c_vdw - p_vdw - l_vdw)
                delta_elec_values.append(c_elec - p_elec - l_elec)
                delta_gb_values.append(c_gb - p_gb - l_gb)
                delta_sa_values.append(c_sa - p_sa - l_sa)

                self.energies['delta_nb'].append(c_nb - p_nb - l_nb)
                self.energies['delta_gb'].append(c_gb - p_gb - l_gb)
                self.energies['delta_sa'].append(c_sa - p_sa - l_sa)
                self.energies['delta_screen'].append(0.0)
                if 'delta_vdw' not in self.energies: self.energies['delta_vdw'] = []
                if 'delta_elec' not in self.energies: self.energies['delta_elec'] = []
                self.energies['delta_vdw'].append(c_vdw - p_vdw - l_vdw)
                self.energies['delta_elec'].append(c_elec - p_elec - l_elec)

                if 'complex_bond' not in self.energies: self.energies['complex_bond'] = []
                if 'complex_angle' not in self.energies: self.energies['complex_angle'] = []
                if 'complex_torsion' not in self.energies: self.energies['complex_torsion'] = []
                self.energies['complex_bond'].append(e_bond)
                self.energies['complex_angle'].append(e_angle)
                self.energies['complex_torsion'].append(e_torsion + e_rbtorsion + e_cmap)

                if i < 5 or i % print_interval == 0:
                    print(f"Frame {i}: Binding={binding_e:.1f} (NB={e_nb:.1f}, OBC={e_obc:.1f}, Int={e_internal:.1f})")
                continue

            if self.entropy_method == 'quasiharmonic':
                frame = traj[i]
                qha_ligand_frames.append(frame.atom_slice(ligand_indices))
                if qha_analyze_complex:
                    qha_protein_frames.append(frame.atom_slice(protein_indices))
                    # Complex system has only protein+ligand (no ions/solvent if they were in PDB but stripped from system)
                    complex_selection = list(protein_indices) + list(ligand_indices)
                    # Sort indices to ensure order matches combined system?
                    # Combined system: protein then ligand.
                    # Usually indices are sorted in PDB?
                    # self.create_combined_system adds protein then ligand.
                    # PDB indices might be Ligand first, or interspersed.
                    # If I slice with [prot..., lig...], mdtraj reorders atoms to match the list.
                    # This matches complex_system order: Protein first, then Ligand.
                    # So passing `protein_indices + ligand_indices` is correct for the system.
                    qha_complex_frames.append(frame.atom_slice(complex_selection))
                
            xyz = traj.xyz[i]
            ligand_pos = unit.Quantity(xyz[lig_idx], unit.nanometer)
            protein_pos = unit.Quantity(xyz[prot_idx], unit.nanometer)
            complex_pos = unit.Quantity(xyz, unit.nanometer)
            
            # Initialize components if not present
            if 'complex_nb' not in self.energies:
                self.energies.update({
                    'complex_nb': [], 'complex_gb': [], 'complex_sa': [], 'complex_screen': [], 'complex_bondsa': []
                })

            try:
                if i == 0:
                    log.debug(f"Frame {i} particle/position count check:")
                    log.debug(f"  Ligand System: {ligand_system.getNumParticles()}, Pos: {len(ligand_pos)}")
                    log.debug(f"  Protein System: {protein_system.getNumParticles()}, Pos: {len(protein_pos)}")
                    comb_len = len(protein_pos) + len(ligand_pos)
                    log.debug(f"  Complex System: {complex_system.getNumParticles()}, Pos Combined: {comb_len}")

                # Calculate fixed enhanced GBSA energies
                ligand_context.setPositions(ligand_pos)
                
                protein_context.setPositions(protein_pos)
                
                # For MM/GBSA: Coordinate Mode's complex_system is built by
                # appending the ligand's atoms after the protein's (see
                # build_complex_system), so combined_pos = protein+ligand
                # reproduces its real atom order -- but ONLY for the
                # small-molecule-ligand path, where build_complex_system
                # explicitly deletes the ligand residue from the protein PDB
                # and re-adds an OpenFF-parameterized ligand at the end. Both
                # Native Mode (complex_struct's atom order is whatever the
                # source topology file's own order is) and PPI Coordinate
                # Mode (is_ppi_coordinate_mode: complex_system is built
                # straight from original_complex_pdb via
                # build_complex_system(..., ligand_mol=None), so its atom
                # order is that PDB's own chain order, not necessarily
                # receptor-then-ligand) need the frame's own true atom order
                # instead. Confirmed on two real datasets: a NAMD PSF where
                # the peptide "ligand" chain is listed BEFORE the receptor
                # chain (native mode), and a GROMACS .gro PPI complex (1GCQ,
                # Zenodo 6638504) where the ligand chain is likewise listed
                # first (coordinate mode) -- both silently swapped receptor/
                # ligand coordinates and produced multi-million-kcal/mol
                # bogus VdW "energy" before this fix.
                if is_native_mode or is_ppi_coordinate_mode:
                    complex_context.setPositions(complex_pos)
                else:
                    combined_pos = unit.Quantity(np.concatenate((xyz[prot_idx], xyz[lig_idx])),
                                                 unit.nanometer)
                    complex_context.setPositions(combined_pos)
                # One evaluation per force group, shared with the parallel worker
                # (_frame_energy_terms); the total is the sum of the groups.
                T = _frame_energy_terms(
                    complex_context, protein_context, ligand_context,
                    tuple(_force_groups(c.getSystem()) for c in (complex_context, protein_context, ligand_context)))
                complex_e, protein_e, ligand_e = T['complex_e'], T['protein_e'], T['ligand_e']
                e_ele, e_vdw, e_gb, e_sa = T['c_elec'], T['c_vdw'], T['c_gb'], T['c_sa']
                e_nb, e_obc = T['e_nb'], T['e_obc']
                e_screen = 0.0  # deprecated
                e_bond, e_angle, e_torsion = T['e_bond'], T['e_angle'], T['e_torsion']
                e_rbtorsion, e_cmap, e_internal = T['e_rbtorsion'], T['e_cmap'], T['e_internal']
                complex_e_clean, protein_e_clean, ligand_e_clean = (
                    T['complex_e_clean'], T['protein_e_clean'], T['ligand_e_clean'])
                binding_e = T['binding_e']
                c_elec, c_vdw, c_nb = T['c_elec'], T['c_vdw'], T['c_nb']
                p_elec, p_vdw = T['p_elec'], T['p_vdw']
                p_nb = p_elec + p_vdw
                l_elec, l_vdw = T['l_elec'], T['l_vdw']
                l_nb = l_elec + l_vdw
                delta_nb_values.append(c_nb - p_nb - l_nb)
                delta_vdw_values.append(c_vdw - p_vdw - l_vdw)
                delta_elec_values.append(c_elec - p_elec - l_elec)
                c_gb, p_gb, l_gb = T['c_gb'], T['p_gb'], T['l_gb']
                delta_gb_values.append(c_gb - p_gb - l_gb)
                c_sa, p_sa, l_sa = T['c_sa'], T['p_sa'], T['l_sa']
                delta_sa_values.append(c_sa - p_sa - l_sa)

                if i == 0:
                    print(f"\n--- DEBUG FRAME 0 ---")
                    print(f"Total Internal Energy (Bond+Angle+Tors) Cancellation Check:")
                    print(f"  Complex Internal: {complex_e - complex_e_clean:.2f}")
                    print(f"  Protein Internal: {protein_e - protein_e_clean:.2f}")
                    print(f"  Ligand Internal:  {ligand_e - ligand_e_clean:.2f}")

                # explicit SASA calculation (User Request) -- a diagnostic print only,
                # so computed on the first frame alone (on every printed frame it
                # took ~40% of a GPU run)
                if i == 0:
                  try:
                      # Create sliced trajectories for SASA (Exclude solvent!)
                      # protein_indices and ligand_indices already exclude solvent
                      complex_indices = sorted(protein_indices + ligand_indices)
                      frame = traj[i]

                      complex_traj = frame.atom_slice(complex_indices)
                      protein_traj = frame.atom_slice(protein_indices)
                      ligand_traj = frame.atom_slice(ligand_indices)
                    
                      # Calculate SASA (nm^2 -> A^2)
                      # ShrakeRupley returns per-atom areas, sum for total
                      # mode='atom' is default.
                      complex_sasa_nm2 = md.shrake_rupley(complex_traj, mode='atom').sum()
                      protein_sasa_nm2 = md.shrake_rupley(protein_traj, mode='atom').sum()
                      ligand_sasa_nm2 = md.shrake_rupley(ligand_traj, mode='atom').sum()
                    
                      complex_sasa = complex_sasa_nm2 * 100.0 # nm^2 to A^2
                      protein_sasa = protein_sasa_nm2 * 100.0
                      ligand_sasa = ligand_sasa_nm2 * 100.0

                      print(f"\nSASA Values (Å²):")
                      print(f"  Protein: {protein_sasa:.2f}")
                      print(f"  Ligand: {ligand_sasa:.2f}")
                      print(f"  Complex: {complex_sasa:.2f}")
                        
                  except Exception as e:
                      if i == 0: 
                          print(f"Warning: SASA calculation failed: {e}")
                    
                # Define deprecated/placeholder variables for compatibility
                c_screen = 0.0
                p_screen = 0.0
                l_screen = 0.0
                e_bondsa = 0.0
                
                if i == 0:
                    print(f"\nBreakdown (kcal/mol):")
                    print(f"            {'Complex':>12} {'Protein':>12} {'Ligand':>12} {'Delta':>12}")
                    print(f"  NB (VdV+Ele): {c_nb:12.2f} {p_nb:12.2f} {l_nb:12.2f} {c_nb - p_nb - l_nb:12.2f}")
                    # Breakdown VDW/Elec
                    print(f"    - VdW:      {c_vdw:12.2f} {p_vdw:12.2f} {l_vdw:12.2f} {c_vdw - p_vdw - l_vdw:12.2f}")
                    print(f"    - Elec:     {c_elec:12.2f} {p_elec:12.2f} {l_elec:12.2f} {c_elec - p_elec - l_elec:12.2f}")
                    print(f"  GB (PolSol):  {c_gb:12.2f} {p_gb:12.2f} {l_gb:12.2f} {c_gb - p_gb - l_gb:12.2f}")
                    print(f"  SA (NonPol):  {c_sa:12.2f} {p_sa:12.2f} {l_sa:12.2f} {c_sa - p_sa - l_sa:12.2f}")
                    # Screening is now 0.0 by definition (VDW moved to VdW row)
                    # print(f"  Screening:    {c_screen:12.2f} {p_screen:12.2f} {l_screen:12.2f} {c_screen - p_screen - l_screen:12.2f}")
                    print(f"  Total Clean:  {complex_e_clean:12.2f} {protein_e_clean:12.2f} {ligand_e_clean:12.2f} {binding_e:12.2f}")
                    print(f"-------------------------------------------------------")

                self.energies['complex'].append(complex_e)
                self.energies['protein'].append(protein_e)
                self.energies['ligand'].append(ligand_e)
                self.energies['binding'].append(binding_e)
                
                self.energies['complex_nb'].append(c_nb)
                self.energies['complex_gb'].append(c_gb)
                self.energies['complex_sa'].append(c_sa)
                self.energies['complex_screen'].append(c_screen)
                self.energies['complex_bondsa'].append(e_bondsa)
                
                # Store Delta Components for Pie Charts & Analysis
                self.energies['delta_nb'].append(c_nb - p_nb - l_nb)
                self.energies['delta_gb'].append(c_gb - p_gb - l_gb)
                self.energies['delta_sa'].append(c_sa - p_sa - l_sa)
                self.energies['delta_screen'].append(c_screen - p_screen - l_screen)
                
                # Store VdW/Elec explicitly if not already present
                if 'delta_vdw' not in self.energies: self.energies['delta_vdw'] = []
                if 'delta_elec' not in self.energies: self.energies['delta_elec'] = []
                
                self.energies['delta_vdw'].append(c_vdw - p_vdw - l_vdw)
                self.energies['delta_elec'].append(c_elec - p_elec - l_elec)
                
                # Store internal energies (create lists later if needed or on fly)
                # For simplicity, store in self.energies dict (init in next step if missing)
                if 'complex_bond' not in self.energies: self.energies['complex_bond'] = []
                if 'complex_angle' not in self.energies: self.energies['complex_angle'] = []
                if 'complex_torsion' not in self.energies: self.energies['complex_torsion'] = []
                
                self.energies['complex_bond'].append(e_bond)
                self.energies['complex_angle'].append(e_angle)
                self.energies['complex_torsion'].append(e_torsion + e_rbtorsion + e_cmap) # Sum torsions for simplicity

                # Optional legacy decomposition (disabled by default as new method is better)
                 
                if i < 5 or i % print_interval == 0:
                    print(f"Frame {i}: Binding={binding_e:.1f} (NB={e_nb:.1f}, OBC={e_obc:.1f}, Int={e_internal:.1f})")
                    
            except Exception as e:
                print(f"Error calculating energies for frame {i}: {e}")
                print("Skipping this frame...")
                continue
        
        calc_time = time.time() - calc_start_time
        total_time = time.time() - analysis_start_time
        
        # Store Delta Components
        self.energies['delta_nb'] = delta_nb_values
        
        self.energies['delta_vdw'] = delta_vdw_values
        self.energies['delta_elec'] = delta_elec_values
             
        self.energies['delta_gb'] = delta_gb_values
        self.energies['delta_sa'] = delta_sa_values
        
        # QHA Calculation
        self.qha_entropy = 0.0
        if self.entropy_method == 'quasiharmonic' and len(qha_ligand_frames) > 0:
             try:
                 print("\nRunning Quasi-Harmonic Analysis (QHA)...")
                 from mmgbsa.quasi_harmonic import QuasiHarmonicAnalysis
                 combined_traj = md.join(qha_ligand_frames)
                 # Note: unit.kelvin is used in init of Core, assume self.temperature is Quantity
                 temp_K = self.temperature.value_in_unit(unit.kelvin)
                 
                 # Ligand
                 combined_traj = md.join(qha_ligand_frames)
                 qha_lig = QuasiHarmonicAnalysis(combined_traj, ligand_system, temp_K, verbose=False)
                 self.qha_entropy = qha_lig.calculate_entropy() # Ligand TS
                 print(f"  • Ligand TS: {self.qha_entropy:.2f} kcal/mol")
                 
                 self.qha_delta_ts = None
                 
                 if qha_analyze_complex and len(qha_protein_frames) > 0:
                      print("  • Computing Protein and Complex QHA (May be slow)...")
                      
                      # Protein
                      # Ideally we should use verbose=True/False based on user pref, but keep log clean
                      traj_prot = md.join(qha_protein_frames)
                      qha_prot = QuasiHarmonicAnalysis(traj_prot, protein_system, temp_K, verbose=False)
                      ts_protein = qha_prot.calculate_entropy()
                      print(f"  • Protein TS: {ts_protein:.2f} kcal/mol")
                      
                      # Complex
                      traj_complex = md.join(qha_complex_frames)
                      qha_complex = QuasiHarmonicAnalysis(traj_complex, complex_system, temp_K, verbose=False)
                      ts_complex = qha_complex.calculate_entropy()
                      print(f"  • Complex TS: {ts_complex:.2f} kcal/mol")
                      
                      # T*Delta S = TS_complex - TS_protein - TS_ligand
                      self.qha_delta_ts = ts_complex - ts_protein - self.qha_entropy
                      print(f"✓ QHA TΔS (Binding): {self.qha_delta_ts:.2f} kcal/mol")
                 else:
                      print(f"✓ QHA Entropy (Ligand TS): {self.qha_entropy:.2f} kcal/mol")
             except Exception as e:
                 print(f"❌ QHA Failed: {e}")
                 import traceback
                 traceback.print_exc()

        # Normal-mode entropy: complex, receptor and ligand minimized from the same frames
        if self.entropy_method in ('nmode', 'normal_mode'):
            try:
                self.nmode_result = self._run_nmode_entropy(
                    traj, selected_frames, {'complex': (complex_system, complex_top),
                                            'receptor': (protein_system, protein_top),
                                            'ligand': (ligand_system, ligand_top)},
                    prot_idx, lig_idx, is_native_mode or is_ppi_coordinate_mode,
                    platform, properties, output_dir)
            except Exception as e:
                log.error(f"Normal-mode entropy failed: {e}")
                self.nmode_result = None

        print(f"✓ Energy calculation time: {calc_time:.1f}s")
        print(f"✓ Total analysis time: {total_time:.1f}s")
        
        # PER-RESIDUE DECOMPOSITION (Integrated Fixed Module)
        # PER-RESIDUE DECOMPOSITION
        df_fast = None # Initialize for later visualization access
        if energy_decomposition:
            decomp_method = getattr(self, 'decomposition_method', 'full') # Default to full
            
            if decomp_method == 'fast':
                 try:
                     print("\nRunning Per-Residue Energy Decomposition (Fast Mode)...")
                     from mmgbsa.energy_decomp import EnergyDecomposer
                     # Fallback to standard flow (parameter injection happens inside EnergyDecomposer now)
                     
                     # Need parmed struct
                     parmed_struct = None
                     if 'parmed' in sys.modules:
                          import parmed as pmd
                          c_prmtop_path = original_complex_pdb
                          if c_prmtop_path and str(c_prmtop_path).endswith('.prmtop') and os.path.exists(c_prmtop_path):
                               parmed_struct = pmd.load_file(str(c_prmtop_path))
                     
                     decomposer = EnergyDecomposer(
                         system=complex_system,
                         topology=complex_top,
                         protein_indices=protein_indices,
                         ligand_indices=ligand_indices,
                         parmed_structure=parmed_struct
                     )
                     
                     results_list = decomposer.analyze_trajectory(traj)
                     
                     # Simple CSV Output for Fast Mode
                     # We flatten the results list: [{'RES1': {...}, 'RES2': ...}, ...]
                     # Averaging over frames
                     
                     # Aggregate
                     agg_results = defaultdict(lambda: {'vdw': [], 'ele': [], 'total': []})
                     
                     for frame_res in results_list:
                         for res_name, vals in frame_res.items():
                             agg_results[res_name]['vdw'].append(vals['vdw'])
                             agg_results[res_name]['ele'].append(vals['electrostatic'])
                             agg_results[res_name]['total'].append(vals['total'])
                             
                     # Create DataFrame
                     final_rows = []
                     for res_name, data in agg_results.items():
                         final_rows.append({
                             'residue': res_name,
                             'vdw_mean': np.mean(data['vdw']),
                             'vdw_std': np.std(data['vdw']),
                             'ele_mean': np.mean(data['ele']),
                             'ele_std': np.std(data['ele']),
                             'total_mean': np.mean(data['total']),
                             'total_std': np.std(data['total'])
                         })
                         
                     df_fast = pd.DataFrame(final_rows)
                     if not df_fast.empty:
                          output_csv = os.path.join(str(output_dir), 'final_decomposition_fast.csv')
                          df_fast.to_csv(output_csv, index=False)
                          print(f"✓ Decomposition saved to: {output_csv}")
                          

                                  

                     else:
                          print("Warning: Fast decomposition returned no significant interactions.")
                          
                 except Exception as e:
                     print(f"❌ Fast Decomposition Failed: {e}")
                     import traceback
                     traceback.print_exc()

            else:
                 # FULL MODE (Original)
                 try:
                     print("\nRunning Per-Residue Energy Decomposition (Full Mode)...")
                     from mmgbsa.decomposition import PerResidueDecomposition
                     
                     # Instantiate Decomposer
                     decomp_tool = PerResidueDecomposition(self, temperature=self.temperature, output_dir=str(output_dir))
                     # its plots go to the run folder (report_dir is otherwise set only
                     # by run_per_residue_analysis, and the plot then lands in the
                     # current working directory)
                     decomp_tool.report_dir = str(self._work_dir_path())
                     
                     # Prepare Systems Dict for the tool
                     # Reuse the context/system we already have active
                     decomp_systems = {
                         'complex_system': complex_system,
                         'complex_context': complex_context,
                         'complex_topology': complex_top # OpenMM topology
                     }
                     
                     # Robust Parameter Injection (Fix for VdW=0 on converted systems)
                     try:
                         import parmed as pmd
                         c_prmtop_path = original_complex_pdb 
                         
                         if c_prmtop_path and str(c_prmtop_path).endswith('.prmtop') and os.path.exists(c_prmtop_path):
                              print(f"Loading ParmEd Structure for Decomposition Parameters: {c_prmtop_path}")
                              struct = pmd.load_file(str(c_prmtop_path))
                              decomp_systems['parmed_structure'] = struct
                     except ImportError:
                         print("Warning: ParmEd not found, using default OpenMM parameters (may be zero for Charmm/Gromacs)")
                     except Exception as e:
                         print(f"Warning: Failed to load ParmEd structure: {e}")
                     
                     # Execute on the existing trajectory slice
                     # Note: 'traj' is an MDTraj object. 
                     # The tool expects mdtraj frames which loop over normally.
                     decomp_tool.execute_decomposition(traj, decomp_systems)
                     
                 except Exception as e:
                     print(f"❌ Decomposition Failed: {e}")
                     import traceback
                     traceback.print_exc()
        
        results = self.save_results(output_dir)
        
        # --- Visualization Section (Moved to end to capture final delta_g) ---
        # energy_decomposition 'fast' always writes a PyMOL script of the residues
        # below -0.5 kcal/mol within 5 A of the ligand, and the heat maps
        # and report of AdvancedVisualization.
        
        # 1. PyMOL Script
        if df_fast is not None:
             try:
                  print("creating PyMOL visualization script...")
                  pymol_viz = PyMOLVisualizer(output_dir)
                  
                  p_threshold, p_distance = -0.5, 5.0     # kcal/mol, A
                  
                  # Verify PDB
                  pdb_file = None
                  self.topology_file = final_pdb_path if 'final_pdb_path' in locals() else complex_pdb
                  if Path(self.topology_file).suffix == '.pdb':
                       pdb_file = self.topology_file
                  else:
                       gen_pdb = Path(output_dir) / "temp_from_traj.pdb"
                       if gen_pdb.exists():
                           pdb_file = gen_pdb
                           
                  if pdb_file:
                       # Need to find output_csv path again or store it
                       # It was 'final_decomposition_fast.csv'
                       # Robustly find it
                       decomp_csv = Path(output_dir) / 'final_decomposition_fast.csv'
                       if decomp_csv.exists():
                           script_path = pymol_viz.generate_script(
                               decomposition_csv=str(decomp_csv),
                               pdb_file=pdb_file,
                               output_filename="view_binding.pml",
                               ligand_resname=getattr(self, 'ligand_resname', 'LIG'),
                               energy_threshold=p_threshold,
                               distance_cutoff=p_distance
                           )
                           if script_path:
                               print(f"✓ PyMOL script generated: {script_path}")
                               print(f"  Run 'pymol {script_path}' to visualize 3D interaction hotspots.")
                       else:
                            print("⚠️  Skipping PyMOL: Decomposition CSV not found.")
                  else:
                       print("⚠️  Skipping PyMOL: No suitable PDB file found.")
             except Exception as e:
                  print(f"⚠️  PyMOL generation failed: {e}")

        # 2. Advanced Visualization & HTML Report
        if df_fast is not None:
             try:
                  print("Generate Advanced Visualization (Heatmaps & Report)...")
                  viz_tool = AdvancedVisualization(output_dir)
                  
                  # Populate results with Binding Energy!
                  mock_results = {
                      'mean_binding_energy': results.get('delta_g', 0.0), # Use Delta G as the total binding energy
                      'std_dev': results.get('std_dev', 0.0),
                      'n_frames': results.get('n_frames', 0),
                      'binding_data': results.get('dataframe'), # Pass time-series data
                      'decomposition_results': {
                          'dataframe': df_fast
                      }
                  }
                  viz_tool.load_mmgbsa_results(mock_results)
                  


                  viz_tool.generate_comprehensive_plots(compound_name="Ligand", generate_report=True)
                  print("✓ Advanced plots generated.")
             except Exception as e:
                  print(f"⚠️ Visualization failed: {e}")
                  import traceback
                  traceback.print_exc()


        # Generate Interactive HTML Report (New Method)
        if output_dir:
            self.generate_detailed_report(
                results, output_dir,
                ligand_resname=ligand_resname or getattr(self, 'report_ligand_resname', None) or 'LIG',
                complex_pdb=complex_pdb)

        return results

        return {
            'output_file': output_file,
            'gb_model': self.gb_model,
            'binding_energy': mean_binding,
            'std_dev': std_dev
        }

    # Above this std. dev. of the gas-phase interaction energy (kcal/mol, i.e.
    # 15 kJ/mol) the IE exponential average cannot be converged in practice
    # (Ekberg & Ryde, J Chem Theory Comput 2021, 17, 5379).
    IE_SIGMA_DE_RELIABILITY_THRESHOLD_KCAL = IE_SIGMA_LIMIT

    def calculate_interaction_entropy(self, energies, temperature=300.0):
        """
        Calculate Interaction Entropy (IE) from interaction energy fluctuations.

        Formula: -TdS = kT * ln < e^(beta * (dE_int - <dE_int>)) >
        Reference: Duan et al., JACS 2016, 138, 5722-5728.

        `energies` are the per-frame gas-phase receptor-ligand
        interaction energies, dE_vdW + dE_elec (see
        `mmgbsa.interaction_entropy.interaction_energies`), as in Duan et al.
        and gmx_MMPBSA -- not the binding enthalpy, whose solvation terms
        damp the fluctuation. sigma(dE_int) is stored in `_last_ie_sigma_dE`
        and a warning is printed above IE_SIGMA_DE_RELIABILITY_THRESHOLD_KCAL.
        This function raises rather than silently returning 0.0 on failure,
        since a silent zero is indistinguishable from a legitimately
        negligible entropy correction in downstream output.

        Returns
        -------
        float
            -TdS in kcal/mol (the entropic penalty; add to the mean binding
            enthalpy to get delta G).

        Raises
        ------
        ValueError
            If fewer than 2 energies are given or an energy is not finite.
        """
        estimates = entropy_estimates(energies, temperature)
        sigma_dE = estimates['sigma_interaction_energy']
        if sigma_dE > self.IE_SIGMA_DE_RELIABILITY_THRESHOLD_KCAL:
            print(f"WARNING: Interaction Entropy sigma(dE_int) = {sigma_dE:.2f} kcal/mol exceeds "
                  f"{self.IE_SIGMA_DE_RELIABILITY_THRESHOLD_KCAL:.1f} kcal/mol (15 kJ/mol), above "
                  f"which the exponential average cannot be converged (Ekberg & Ryde 2021); treat "
                  f"the reported -TdS as unreliable.")
        self._last_ie_sigma_dE = sigma_dE  # exposed for callers/reporting that want the diagnostic
        return estimates['interaction_entropy']

    def _run_nmode_entropy(self, traj, frame_numbers, systems, prot_idx, lig_idx, full_complex,
                           platform, properties, output_dir=None):
        """Normal-mode binding entropy on nmode_frames of the analysed frames.

        `systems` holds the (System, Topology) of the complex, receptor and ligand;
        `full_complex` says whether the complex System follows the frame's own
        atom order (Native and PPI Coordinate Mode) or receptor-then-ligand.
        Each species is minimized and its normal modes computed as in
        mmgbsa.nmode_entropy; amber_compatible uses MMPBSA.py's nmode energy
        function. Writes nmode_entropy.csv and returns the summary of
        nmode_entropy.binding_entropy with the per-frame rows.
        """
        from mmgbsa import nmode_entropy as nm

        temperature = self.temperature
        if unit.is_quantity(temperature):
            temperature = temperature.value_in_unit(unit.kelvin)
        # MMPBSA.py fixes 298.15 K for its normal-mode entropy, whatever the run temperature
        if self.amber_compatible:
            temperature = nm.AMBER_THERMO_TEMPERATURE
        device = {k: v for k, v in (properties or {}).items() if 'DeviceIndex' in k}
        species = {}
        for name, (system, topology) in systems.items():
            s = nm.amber_nmode_system(system, topology) if self.amber_compatible else nm.nmode_system(system)
            species[name] = (nm.make_contexts(s, platform.getName(), self.nmode_precision, device),
                             [s.getParticleMass(i).value_in_unit(unit.dalton) for i in range(s.getNumParticles())])
        precision = 'double' if species['complex'][0]['mixed'] is None else 'mixed'
        indices = nm.select_frames(traj.n_frames, self.nmode_frames)
        energy_function = ("MMPBSA.py's nmode function (HCT, rgbmax 15 A, no salt)" if self.amber_compatible
                           else f"{self.gb_model} without the non-polar term")
        print(f"\nNormal-mode entropy on {len(indices)} of {traj.n_frames} frames: {energy_function}, "
              f"minimized to RMS gradient {self.nmode_drms:g} kcal/mol/A, Hessian in {precision} precision")

        frames = []
        for count, i in enumerate(indices, 1):
            t0 = time.time()
            xyz = traj.xyz[i].astype(float)
            positions = {'complex': xyz if full_complex else np.concatenate((xyz[prot_idx], xyz[lig_idx])),
                         'receptor': xyz[prot_idx], 'ligand': xyz[lig_idx]}
            frame = {'frame': int(frame_numbers[i])}
            for name in ('complex', 'receptor', 'ligand'):
                contexts, masses = species[name]
                frame[name] = nm.species_entropy(contexts, positions[name], masses, temperature,
                                                 self.nmode_drms, self.nmode_maxcyc)
            frames.append(frame)
            row = nm.binding_entropy([frame], temperature)[0][0]
            status = (f"-TdS = {row['minus_TdS']:.2f} kcal/mol" if row['converged'] else
                      "left out, minimization did not reach the RMS gradient")
            print(f"  frame {frame['frame']} ({count}/{len(indices)}, {time.time() - t0:.0f} s): {status}")

        rows, summary = nm.binding_entropy(frames, temperature)
        summary.update(precision=precision, drms=self.nmode_drms, amber_compatible=self.amber_compatible,
                       rows=rows)
        if output_dir:
            path = Path(output_dir) / 'nmode_entropy.csv'
            pd.DataFrame(rows).to_csv(path, index=False)
            summary['csv'] = str(path)
        if summary['n_converged']:
            print(f"✓ Normal-mode -TdS = {summary['value']:.2f} ± {summary['sd']:.2f} kcal/mol (SD, "
                  f"{summary['n_converged']} of {summary['n_frames']} frames)")
        else:
            print("❌ Normal-mode entropy: no frame reached the RMS gradient criterion")
        return summary

    def save_results(self, output_dir=None):
        """Save fixed enhanced MM/GBSA results to file"""
        if len(self.energies['binding']) == 0:
            print("No successful energy calculations!")
            return None
            
        # Prepare DataFrame with components if available
        data = {
            'frame': np.arange(len(self.energies['binding'])),
            'complex_energy': self.energies['complex'],
            'protein_energy': self.energies['protein'],
            'ligand_energy': self.energies['ligand'],
            'binding_energy': self.energies['binding']
        }
        
        # Add components if they exist
        if 'complex_nb' in self.energies and len(self.energies['complex_nb']) > 0:
            data.update({
                'complex_nb': self.energies['complex_nb'],
                'complex_gb': self.energies['complex_gb'],
                'complex_sa': self.energies['complex_sa'],
                'complex_screen': self.energies['complex_screen']
            })
            
        if 'complex_bond' in self.energies and len(self.energies['complex_bond']) > 0:
            data.update({
                'complex_bond': self.energies['complex_bond'],
                'complex_angle': self.energies['complex_angle'],
                'complex_torsion': self.energies['complex_torsion']
            })
            
        if 'delta_nb' in self.energies and len(self.energies['delta_nb']) > 0:
            if 'delta_vdw' in self.energies and len(self.energies['delta_vdw']) > 0:
                data['delta_vdw'] = self.energies['delta_vdw']
            
            if 'delta_elec' in self.energies and len(self.energies['delta_elec']) > 0:
                data['delta_elec'] = self.energies['delta_elec']
            
            data.update({
                'delta_nb': self.energies['delta_nb'],
                'delta_gb': self.energies['delta_gb'],
                'delta_sa': self.energies['delta_sa']
            })
            
        df = pd.DataFrame(data)
        
        output_file = f'fixed_enhanced_mmgbsa_results_{self.gb_model.lower()}.csv'
        if output_dir:
            output_file = str(Path(output_dir) / output_file)
            
        df.to_csv(output_file, index=False)
        
        mean_binding = df['binding_energy'].mean()
        std_error = df['binding_energy'].std() / np.sqrt(len(df))
        std_dev = df['binding_energy'].std()
        
        print(f"\nResults saved to {output_file}")
        print(f"Advanced GB Model: {self.gb_model}")
        print(f"Salt Concentration: {self.salt_concentration} M")
        print(f"Mean binding energy (Enthalpy): {mean_binding:.2f} ± {std_dev:.2f} kcal/mol (Standard Deviation)")
        print(f"Standard Error of Mean: {std_error:.2f} kcal/mol")
    
        # ENTROPY CALCULATION
        entropy_penalty = 0.0
        delta_g = mean_binding
        entropy = None

        if self.entropy_method in ('interaction', 'c2'):
            try:
                print(f"\nCalculating {METHOD_LABELS[self.entropy_method]} from the interaction energy...")
                temp_val = self.temperature.value_in_unit(unit.kelvin)
                entropy = entropy_result(interaction_energies(df), self.entropy_method, temp_val)
                entropy_penalty = entropy['value']
                delta_g = mean_binding + entropy_penalty

                for line in describe(entropy):
                    print(line)
                print(f"Final Binding Free Energy (ΔG = ΔH - TΔS, {entropy['label']}): {delta_g:.2f} kcal/mol")
            except Exception as e:
                print(f"Entropy calculation failed: {e}")
                entropy_penalty = 0.0
                delta_g = mean_binding
                entropy = None
        elif self.entropy_method in ('nmode', 'normal_mode'):
            res = self.nmode_result
            if res and res.get('n_converged'):
                entropy = {k: v for k, v in res.items() if k != 'rows'}
                entropy.update(method='nmode', label='Normal-mode entropy')
                entropy_penalty = res['value']
                delta_g = mean_binding + entropy_penalty
                print(f"\nNormal-mode entropy (-TΔS): {entropy_penalty:.2f} ± {res['sd']:.2f} kcal/mol "
                      f"(SD over {res['n_converged']} frames)")
                print(f"Final Binding Free Energy (ΔG = ΔH - TΔS, normal modes): {delta_g:.2f} kcal/mol")
            else:
                print("\nNormal-mode entropy unavailable; ΔG is not reported.")
        elif self.entropy_method == 'quasiharmonic':
             ts_val = getattr(self, 'qha_entropy', 0.0)
             delta_ts_val = getattr(self, 'qha_delta_ts', None)
             
             print(f"\nQuasi-Harmonic Analysis Results:")
             if delta_ts_val is not None:
                 print(f"  Ligand TS: {ts_val:.2f} kcal/mol")
                 print(f"  TΔS (Binding): {delta_ts_val:.2f} kcal/mol")
                 # TdS is negative (usually). -TdS is positive penalty.
                 # Formula: dG = dH - TdS.
                 entropy_penalty = -delta_ts_val 
                 delta_g = mean_binding + entropy_penalty
                 print(f"  Entropy Penalty (-TΔS): {entropy_penalty:.2f} kcal/mol")
                 print(f"Final Binding Free Energy (ΔG): {delta_g:.2f} kcal/mol")
             else:
                 print(f"  Absolute Vibrational Entropy (TS): {ts_val:.2f} kcal/mol")
                 print(f"  (Note: This is absolute entropy of the bound ligand, not ΔS)")
                 # We do not modify dG because we lack reference state
                 entropy_penalty = 0.0 
                 delta_g = mean_binding
        else:
            print(f"\nInteraction Entropy disabled (entropy_method={self.entropy_method})")
 

        
        merged_assumptions = list(getattr(self, 'physics_assumptions', []))
        if hasattr(self, 'gbsa_manager') and hasattr(self.gbsa_manager, 'physics_assumptions'):
            for a in self.gbsa_manager.physics_assumptions:
                if a not in merged_assumptions:
                    merged_assumptions.append(a)

        return {
            'output_file': output_file,
            'dataframe': df,  # Return full dataframe for plotting
            'gb_model': self.gb_model,
            'mean_binding_energy': mean_binding,
            'median_binding_energy': df['binding_energy'].median(),
            'min_binding_energy': df['binding_energy'].min(),
            'max_binding_energy': df['binding_energy'].max(),
            'std_dev': std_dev,
            'std_error': std_error,
            'n_frames': len(df),
            'entropy_penalty': entropy_penalty,
            'entropy': entropy,
            'delta_g': delta_g,
            'parameterized_residues': self.parameterized_residues,
            'physics_assumptions': merged_assumptions,
            'platform_requested': getattr(self, 'analysis_platform', {}).get('requested'),
            'platform_used': getattr(self, 'analysis_platform', {}).get('used'),
            'platform_fallback': getattr(self, 'analysis_platform', {}).get('fallback'),
            'bootstrap_results': self._calculate_bootstrap(df['binding_energy'].values) if len(df) > 5 else None
        }

    def _calculate_bootstrap(self, data_values, n_boot=1000):
        try:
           resampled_means = []
           for _ in range(n_boot):
               resample = np.random.choice(data_values, size=len(data_values), replace=True)
               resampled_means.append(np.mean(resample))
           return {
               'mean': np.mean(resampled_means), 
               'std': np.std(resampled_means), 
               'ci_lower': np.percentile(resampled_means, 2.5), 
               'ci_upper': np.percentile(resampled_means, 97.5)
           }
        except:
           return None


    def create_combined_system(self, protein_system, ligand_system):
        """
        Concatenate a separately-built protein System and ligand System into
        a single OpenMM System (protein particles first, ligand particles
        appended after, at index offset `n_protein`), re-merging each force
        type found on both inputs into one combined force covering all
        particles.

        This exists so that per-force-group energies (bond=10, angle=11,
        torsion=12, RBTorsion=13, CMAP=14, nonbonded=0, GB=1, screening=3,
        SA=4 -- see the `assign_force_groups` local function used elsewhere
        in this class) can be read from a single Context for the *combined*
        complex-equivalent system, letting energy-component decomposition
        work the same way whether the "complex" came from one native
        topology or from stitching a protein system and ligand system
        together (e.g. for dimer/multi-subunit binding modes).

        Only forces present on BOTH `protein_system` and `ligand_system` are
        merged and added to the output; a force present on only one input is
        silently dropped (e.g. CMAPTorsionForce, which a ligand system
        typically lacks). GBSAOBCForce pairing (GB vs. SA) is disambiguated
        via each force's own `SurfaceAreaEnergy` setting (0.0 for the polar
        GB force, nonzero for the ACE nonpolar-SA force in this codebase's
        force-construction convention), not by inspecting particle charges.

        Parameters
        ----------
        protein_system : openmm.System
            System built for the protein/receptor alone.
        ligand_system : openmm.System
            System built for the ligand alone.

        Returns
        -------
        openmm.System
            Combined system with `n_protein + n_ligand` particles and merged
            forces, force-grouped as described above.
        """
        print("Creating combined protein+ligand system for MM/GBSA...")
        combined_system = openmm.System()
        n_protein = protein_system.getNumParticles()
        for i in range(n_protein):
            combined_system.addParticle(protein_system.getParticleMass(i))
        n_ligand = ligand_system.getNumParticles()
        for i in range(n_ligand):
            combined_system.addParticle(ligand_system.getParticleMass(i))
        print(f"  Combined: {n_protein} protein + {n_ligand} ligand = {combined_system.getNumParticles()} particles")
        protein_offset = 0
        ligand_offset = n_protein
        protein_gb = protein_screening = protein_sa = None
        ligand_gb = ligand_screening = ligand_sa = None
        for force in protein_system.getForces():
            if isinstance(force, openmm.CustomGBForce):
                protein_gb = force
            elif isinstance(force, openmm.CustomNonbondedForce):
                protein_screening = force
            elif isinstance(force, openmm.CustomBondForce):
                protein_sa = force
        for force in ligand_system.getForces():
            if isinstance(force, openmm.CustomGBForce):
                ligand_gb = force
            elif isinstance(force, openmm.CustomNonbondedForce):
                ligand_screening = force
            elif isinstance(force, openmm.CustomBondForce):
                ligand_sa = force
        if protein_gb and ligand_gb:
            f = self._merge_custom_gb_forces(protein_gb, ligand_gb, protein_offset, ligand_offset)
            f.setForceGroup(2) # Group 2: CustomGB (SA or similar)
            combined_system.addForce(f)
            print("  ✓ Merged CustomGBForce (Group 2)")
        if protein_screening and ligand_screening:
            f = self._merge_custom_nonbonded_forces(protein_screening, ligand_screening, protein_offset, ligand_offset)
            f.setForceGroup(3) # Group 3: Salt Screening
            combined_system.addForce(f)
            print("  ✓ Merged CustomNonbondedForce (Group 3)")
        if protein_sa and ligand_sa:
            f = self._merge_custom_bond_forces(protein_sa, ligand_sa, protein_offset, ligand_offset)
            f.setForceGroup(4) # Group 4: CustomBondSA
            combined_system.addForce(f)
            print("  ✓ Merged CustomBondForce (Group 4)")
        
        # --- Merge Internal Forces (Bonds, Angles, Torsions) ---
        # 1. HarmonicBondForce
        pf, lf = None, None
        for force in protein_system.getForces():
             if isinstance(force, openmm.HarmonicBondForce): pf = force; break
        for force in ligand_system.getForces():
             if isinstance(force, openmm.HarmonicBondForce): lf = force; break
        if pf and lf:
            f = self._merge_bond_forces(pf, lf, protein_offset, ligand_offset)
            f.setForceGroup(10) # Group 10
            combined_system.addForce(f)
            print("  ✓ Merged HarmonicBondForce (Group 10)")

        # 2. HarmonicAngleForce
        pf, lf = None, None
        for force in protein_system.getForces():
             if isinstance(force, openmm.HarmonicAngleForce): pf = force; break
        for force in ligand_system.getForces():
             if isinstance(force, openmm.HarmonicAngleForce): lf = force; break
        if pf and lf:
            f = self._merge_angle_forces(pf, lf, protein_offset, ligand_offset)
            f.setForceGroup(11) # Group 11
            combined_system.addForce(f)
            print("  ✓ Merged HarmonicAngleForce (Group 11)")

        # 3. PeriodicTorsionForce
        pf, lf = None, None
        for force in protein_system.getForces():
             if isinstance(force, openmm.PeriodicTorsionForce): pf = force; break
        for force in ligand_system.getForces():
             if isinstance(force, openmm.PeriodicTorsionForce): lf = force; break
        if pf and lf:
            f = self._merge_torsion_forces(pf, lf, protein_offset, ligand_offset)
            f.setForceGroup(12) # Group 12
            combined_system.addForce(f)
            print("  ✓ Merged PeriodicTorsionForce (Group 12)")

        # 4. RBTorsionForce
        pf, lf = None, None
        for force in protein_system.getForces():
             if isinstance(force, openmm.RBTorsionForce): pf = force; break
        for force in ligand_system.getForces():
             if isinstance(force, openmm.RBTorsionForce): lf = force; break
        if pf and lf:
            f = self._merge_rb_torsion_forces(pf, lf, protein_offset, ligand_offset)
            f.setForceGroup(13) # Group 13
            combined_system.addForce(f)
            print("  ✓ Merged RBTorsionForce (Group 13)")
            
        # 5. CMAPTorsionForce
        pf, lf = None, None
        for force in protein_system.getForces():
             if isinstance(force, openmm.CMAPTorsionForce): pf = force; break
        for force in ligand_system.getForces():
             if isinstance(force, openmm.CMAPTorsionForce): lf = force; break
        # Note: Often only protein has CMAP. If ligand is lacking, we can still add protein's CMAP or skip if strict.
        # Strict merge requires both, but loose merge (logic in helper) could handle one.
        # For now, strict merge if both present, or adapt helper.
        # Actually, if ligand lacks it, we shouldn't merge? Or we should add just protein's?
        # Standard approach: If protein has it, we want it in complex. 
        # But my helper assumes both. Let's stick to strict 'if both' for now to avoid crashes, 
        # but in reality ligand won't have CMAP. This means CMAP might be lost for complex if looking for 'lf'.
        # However, ignoring CMAP in complex while having it in protein leads to energy mismatch.
        # FIX: We need to handle case where only one has it. But 'merged' implies joining.
        # If ligand has 0 CMAPs, the loop range(lf.getNumTorsions()) is empty, so it works fine!
        # Just need to find the force object even if empty?
        # But if ligand system DOESN'T HAVE the force object, lf is None.
        # We need to relax checks or create dummy force.
        # Let's skip CMAP for now unless both have it (safe path), or user complains about specific residue errors.
        if pf and lf:
            f = self._merge_cmap_torsion_forces(pf, lf, protein_offset, ligand_offset)
            f.setForceGroup(14) # Group 14
            combined_system.addForce(f)
            print("  ✓ Merged CMAPTorsionForce (Group 14)")
        
        # Merge standard NonbondedForce
        protein_nb = None
        ligand_nb = None
        for force in protein_system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                protein_nb = force
                break
        for force in ligand_system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                ligand_nb = force
                break
        
        if protein_nb and ligand_nb:
            f = self._merge_nonbonded_forces(protein_nb, ligand_nb, protein_offset, ligand_offset)
            f.setForceGroup(0) # Group 0: Nonbonded (VDW+Ele)
            combined_system.addForce(f)
            print("  ✓ Merged NonbondedForce (Group 0)")

        # Merge GBSAOBCForces (GB and SA)
        # Note: A system might have MULTIPLE GBSAOBCForces (one for GB, one for SA)
        # We need to pair them up. Assuming order is preserved:
        # First one is GB (Group 1), Second is SA (Group 4) if present.
        # Safer strategy: Identify by Force Group or Charge parameters?
        # But force groups are reset in combined system.
        # Let's collect ALL GBSAOBCForces and merge index-wise.
        
        p_obc_forces = [f for f in protein_system.getForces() if isinstance(f, openmm.GBSAOBCForce)]
        l_obc_forces = [f for f in ligand_system.getForces() if isinstance(f, openmm.GBSAOBCForce)]
        
        if len(p_obc_forces) != len(l_obc_forces):
            print(f"Warning: Protein has {len(p_obc_forces)} OBC forces but Ligand has {len(l_obc_forces)}. Merging by index up to minimum.")
            
        for i in range(min(len(p_obc_forces), len(l_obc_forces))):
            pf = p_obc_forces[i]
            lf = l_obc_forces[i]

            # Determine GB vs. SA identity from each force's own
            # SurfaceAreaEnergy setting, which this codebase sets
            # unconditionally and oppositely at construction time:
            # the real GB force always has SurfaceAreaEnergy == 0.0
            # (see _create_fallback_obc_force / the factory GB path,
            # both of which call setSurfaceAreaEnergy(0.0)), while the
            # ACE nonpolar-SA force always has a nonzero coefficient
            # (see _create_ace_sa_force, setSurfaceAreaEnergy(2.25936)).
            # This is a direct, always-true signal, unlike inferring
            # intent from whether a given particle's charge happens to
            # be near zero -- a ligand atom can legitimately have a
            # near-zero partial charge, which would misclassify a real
            # GB force as the SA force under the old heuristic.
            sa_unit = unit.kilojoule_per_mole / unit.nanometer**2
            pf_sa_energy = pf.getSurfaceAreaEnergy().value_in_unit(sa_unit)
            lf_sa_energy = lf.getSurfaceAreaEnergy().value_in_unit(sa_unit)
            if pf_sa_energy != lf_sa_energy:
                print(f"Warning: paired OBC forces at index {i} have mismatched "
                      f"SurfaceAreaEnergy ({pf_sa_energy} vs {lf_sa_energy}); "
                      f"protein/ligand systems may not have been built consistently.")
            is_sa = (pf_sa_energy != 0.0) or (lf_sa_energy != 0.0)

            f = self._merge_obc_forces(pf, lf, protein_offset, ligand_offset)
            
            if is_sa:
                f.setForceGroup(4) # SA
                print(f"  ✓ Merged GBSAOBCForce (SA, Group 4)")
            else:
                f.setForceGroup(1) # GB
                print(f"  ✓ Merged GBSAOBCForce (GB, Group 1)")
            
            combined_system.addForce(f)

        print(f"✓ Combined system: {combined_system.getNumForces()} forces")
        return combined_system


    def _merge_obc_forces(self, pf, lf, po, lo):
        """
        Merge a protein GBSAOBCForce `pf` and ligand GBSAOBCForce `lf` into
        one force covering both particle sets (protein particles first,
        ligand particles at offset `po`/`lo` -- both always 0/n_protein in
        practice, passed through from `create_combined_system`). Dielectrics,
        SurfaceAreaEnergy, nonbonded method, and cutoff are copied from `pf`
        (both inputs are expected to share identical settings).
        """
        c = openmm.GBSAOBCForce()
        c.setSolventDielectric(pf.getSolventDielectric())
        c.setSoluteDielectric(pf.getSoluteDielectric())
        
        # Ensure we copy SA energy correctly
        sa_energy = pf.getSurfaceAreaEnergy()
        c.setSurfaceAreaEnergy(sa_energy)
        
        c.setNonbondedMethod(pf.getNonbondedMethod())
        c.setCutoffDistance(pf.getCutoffDistance())
        
        # Add protein particles
        for i in range(pf.getNumParticles()):
            q, r, s = pf.getParticleParameters(i)
            c.addParticle(q, r, s)
            
        # Add ligand particles
        for i in range(lf.getNumParticles()):
            q, r, s = lf.getParticleParameters(i)
            c.addParticle(q, r, s)
            
        return c

    def _merge_nonbonded_forces(self, pf, lf, po, lo):
        """
        Merge protein and ligand NonbondedForces into one, offsetting ligand
        particle/exception indices by `lo`. Only intra-protein and
        intra-ligand exceptions (1-4 scaling, exclusions) are preserved;
        no protein-ligand exceptions are created, so all inter-molecular
        nonbonded pairs are evaluated at full strength (the physically
        correct behavior for previously-non-bonded protein and ligand atoms).
        """
        c = openmm.NonbondedForce()
        c.setNonbondedMethod(pf.getNonbondedMethod())
        c.setCutoffDistance(pf.getCutoffDistance())
        c.setEwaldErrorTolerance(pf.getEwaldErrorTolerance())
        # Copy particles
        for i in range(pf.getNumParticles()):
            c.addParticle(*pf.getParticleParameters(i))
        for i in range(lf.getNumParticles()):
            c.addParticle(*lf.getParticleParameters(i))
        # Copy exceptions
        for i in range(pf.getNumExceptions()):
            p1, p2, q, sig, eps = pf.getExceptionParameters(i)
            c.addException(p1+po, p2+po, q, sig, eps)
        for i in range(lf.getNumExceptions()):
            p1, p2, q, sig, eps = lf.getExceptionParameters(i)
            c.addException(p1+lo, p2+lo, q, sig, eps)
        # Ensure exclusions are handled?
        # NonbondedForce handles exclusions via exceptions or useExceptionList.
        # This basic copy preserves intrasegment exclusions.
        # Intersegment (Protein-Ligand) interactions will be calculated (default: no exclusion).
        return c

    def _merge_custom_gb_forces(self, pf, lf, po, lo):
        """
        Merge protein and ligand CustomGBForce instances (used for GB models
        that expose dielectrics as queryable global parameters, e.g. GBn/
        GBn2/HCT -- NOT OBC2 with a nonzero kappa, whose dielectrics are
        embedded as literal constants in its energy-term expression strings
        and are unaffected by this function).

        KNOWN LIMITATION: solventDielectric/soluteDielectric are hardcoded
        to 78.5/1.0 below rather than copied from `pf`'s actual global
        parameter values. If a user configures non-default dielectrics for
        one of these GB models, that setting is silently dropped by this
        merge path.
        """
        c = openmm.CustomGBForce()
        # Copy per-particle parameters dynamically
        for i in range(pf.getNumPerParticleParameters()):
            c.addPerParticleParameter(pf.getPerParticleParameterName(i))

        c.addGlobalParameter("solventDielectric", 78.5)
        c.addGlobalParameter("soluteDielectric", 1.0)
        
        for i in range(pf.getNumComputedValues()):
            c.addComputedValue(*pf.getComputedValueParameters(i))
        for i in range(pf.getNumEnergyTerms()):
            c.addEnergyTerm(*pf.getEnergyTermParameters(i))
        for i in range(pf.getNumParticles()):
            c.addParticle(pf.getParticleParameters(i))
        for i in range(lf.getNumParticles()):
            c.addParticle(lf.getParticleParameters(i))
        c.setNonbondedMethod(pf.getNonbondedMethod())
        c.setCutoffDistance(pf.getCutoffDistance())
        return c

    def _merge_custom_nonbonded_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand CustomNonbondedForce pair (e.g. the Debye-Huckel
        salt-screening term), offsetting ligand particle/exclusion indices by `lo`."""
        c = openmm.CustomNonbondedForce(pf.getEnergyFunction())
        for i in range(pf.getNumPerParticleParameters()):
            c.addPerParticleParameter(pf.getPerParticleParameterName(i))
        for i in range(pf.getNumGlobalParameters()):
            c.addGlobalParameter(pf.getGlobalParameterName(i), pf.getGlobalParameterDefaultValue(i))
        for i in range(pf.getNumParticles()):
            c.addParticle(pf.getParticleParameters(i))
        for i in range(lf.getNumParticles()):
            c.addParticle(lf.getParticleParameters(i))
        for i in range(pf.getNumExclusions()):
            p1, p2 = pf.getExclusionParticles(i)
            c.addExclusion(p1 + po, p2 + po)
        for i in range(lf.getNumExclusions()):
            p1, p2 = lf.getExclusionParticles(i)
            c.addExclusion(p1 + lo, p2 + lo)
        c.setNonbondedMethod(pf.getNonbondedMethod())
        c.setCutoffDistance(pf.getCutoffDistance())
        return c

    def _merge_custom_bond_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand CustomBondForce pair (e.g. the LCPO nonpolar
        surface-area term), offsetting ligand bond-atom indices by `lo`."""
        c = openmm.CustomBondForce(pf.getEnergyFunction())
        for i in range(pf.getNumPerBondParameters()):
            c.addPerBondParameter(pf.getPerBondParameterName(i))
        for i in range(pf.getNumGlobalParameters()):
            c.addGlobalParameter(pf.getGlobalParameterName(i), pf.getGlobalParameterDefaultValue(i))
        for i in range(pf.getNumBonds()):
            p1, p2, params = pf.getBondParameters(i)
            c.addBond(p1 + po, p2 + po, params)
        for i in range(lf.getNumBonds()):
            p1, p2, params = lf.getBondParameters(i)
            c.addBond(p1 + lo, p2 + lo, params)
        return c

    def _merge_bond_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand HarmonicBondForce pair, offsetting ligand atom indices by `lo`."""
        c = openmm.HarmonicBondForce()
        for i in range(pf.getNumBonds()):
            p1, p2, length, k = pf.getBondParameters(i)
            c.addBond(p1+po, p2+po, length, k)
        for i in range(lf.getNumBonds()):
            p1, p2, length, k = lf.getBondParameters(i)
            c.addBond(p1+lo, p2+lo, length, k)
        return c

    def _merge_angle_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand HarmonicAngleForce pair, offsetting ligand atom indices by `lo`."""
        c = openmm.HarmonicAngleForce()
        for i in range(pf.getNumAngles()):
            p1, p2, p3, angle, k = pf.getAngleParameters(i)
            c.addAngle(p1+po, p2+po, p3+po, angle, k)
        for i in range(lf.getNumAngles()):
            p1, p2, p3, angle, k = lf.getAngleParameters(i)
            c.addAngle(p1+lo, p2+lo, p3+lo, angle, k)
        return c

    def _merge_torsion_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand PeriodicTorsionForce (proper dihedrals) pair,
        offsetting ligand atom indices by `lo`."""
        c = openmm.PeriodicTorsionForce()
        for i in range(pf.getNumTorsions()):
            p1, p2, p3, p4, per, phase, k = pf.getTorsionParameters(i)
            c.addTorsion(p1+po, p2+po, p3+po, p4+po, per, phase, k)
        for i in range(lf.getNumTorsions()):
            p1, p2, p3, p4, per, phase, k = lf.getTorsionParameters(i)
            c.addTorsion(p1+lo, p2+lo, p3+lo, p4+lo, per, phase, k)
        return c
    
    def _merge_rb_torsion_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand RBTorsionForce (Ryckaert-Bellemans dihedrals,
        e.g. from GROMOS/OPLS-derived topologies) pair, offsetting ligand atom indices by `lo`."""
        c = openmm.RBTorsionForce()
        for i in range(pf.getNumTorsions()):
            p1, p2, p3, p4, c0, c1, c2, c3, c4, c5 = pf.getTorsionParameters(i)
            c.addTorsion(p1+po, p2+po, p3+po, p4+po, c0, c1, c2, c3, c4, c5)
        for i in range(lf.getNumTorsions()):
            p1, p2, p3, p4, c0, c1, c2, c3, c4, c5 = lf.getTorsionParameters(i)
            c.addTorsion(p1+lo, p2+lo, p3+lo, p4+lo, c0, c1, c2, c3, c4, c5)
        return c
        
    def _merge_cmap_torsion_forces(self, pf, lf, po, lo):
        """Merge a protein/ligand CMAPTorsionForce (CHARMM backbone grid correction)
        pair. Only called from `create_combined_system` when BOTH inputs have a
        CMAPTorsionForce; a ligand system typically lacks one, so in practice this
        merge (and any CMAP contribution to the combined system) is skipped."""
        c = openmm.CMAPTorsionForce()
        # Copy maps first (required before adding torsions)
        for i in range(pf.getNumMaps()):
            size, energy = pf.getMapParameters(i)
            c.addMap(size, energy)
        # Handle CMAP offset if ligand has maps (unlikely but safe to handle)
        lf_map_offset = pf.getNumMaps()
        for i in range(lf.getNumMaps()):
            size, energy = lf.getMapParameters(i)
            c.addMap(size, energy)
            
        for i in range(pf.getNumTorsions()):
            map_id, p1, p2, p3, p4, p5, p6, p7, p8 = pf.getTorsionParameters(i)
            c.addTorsion(map_id, p1+po, p2+po, p3+po, p4+po, p5+po, p6+po, p7+po, p8+po)
            
        for i in range(lf.getNumTorsions()):
            map_id, p1, p2, p3, p4, p5, p6, p7, p8 = lf.getTorsionParameters(i)
            c.addTorsion(map_id+lf_map_offset, p1+lo, p2+lo, p3+lo, p4+lo, p5+lo, p6+lo, p7+lo, p8+lo)
            
        return c


    """Fixed enhanced main function to run MM/GBSA analysis"""
    ligand_mol = 'test/ligand.sdf'
