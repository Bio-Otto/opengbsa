"""
Ligand and protein parameterization helpers.
Provides a ParameterizationManager that wraps common parameterization tasks.
"""

import os
import warnings

try:
    from openff.toolkit.topology import Molecule
    OPENFF_AVAILABLE = True
except Exception:
    OPENFF_AVAILABLE = False

try:
    import parmed as pmd
    PARMED_AVAILABLE = True
except Exception:
    PARMED_AVAILABLE = False


class ParameterizationManager:
    """Manage parameterization of ligands and proteins."""

    def __init__(self, verbose=False):
        self.verbose = verbose

    def parameterize_ligand_openff(self, ligand_file, charge_method='am1bcc'):
        """
        Parameterize ligand using OpenFF if available.
        Returns a ParmEd Structure-like object or None.
        """
        if not OPENFF_AVAILABLE:
            warnings.warn("OpenFF toolkit not available; ligand parameterization skipped.")
            return None

        try:
            mol = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
            if self.verbose:
                print(f"✓ Loaded ligand via OpenFF: {ligand_file}")

            # Charge assignment (best-effort; real implementation may vary)
            if charge_method and hasattr(mol, 'assign_partial_charges'):
                try:
                    mol.assign_partial_charges(charge_method)
                except Exception:
                    pass

            # Convert to ParmEd Structure if ParmEd available
            if PARMED_AVAILABLE:
                # Use temporary conversion via mol.to_file if needed; keep simple here
                tmp_pdb = ligand_file + '.openff.tmp.pdb'
                mol.to_file(tmp_pdb)
                struct = pmd.load_file(tmp_pdb)
                try:
                    os.remove(tmp_pdb)
                except Exception:
                    pass
                return struct
            else:
                return mol
        except Exception as e:
            if self.verbose:
                print(f"⚠ Ligand parameterization failed: {e}")
            return None

    def parameterize_protein_amber(self, protein_pdb, forcefield='amber99sb-ildn'):
        """
        Parameterize protein using OpenMM ForceField (best-effort wrapper).
        Returns an OpenMM-compatible ForceField or ParmEd Structure.
        """
        try:
            import openmm.app as app
            ff = app.ForceField(f'{forcefield}.xml')
            if self.verbose:
                print(f"✓ Loaded forcefield: {forcefield}")
            return ff
        except Exception as e:
            if self.verbose:
                print(f"⚠ Protein parameterization failed: {e}")
            return None

    def set_ligand_forcefield(self, top, ligand_struct, ff_engine='openmm'):
        """
        Attach ligand parameters to a topology. Minimal placeholder.
        """
        if self.verbose:
            print("• set_ligand_forcefield called (placeholder)")
        return True

    def get_ligand_positions_openmm(self, ligand_struct):
        """
        Return ligand positions in OpenMM-friendly format (list of Vec3).
        Accepts either a ParmEd Structure or an OpenFF Molecule-like object.
        """
        try:
            from simtk import unit
        except Exception:
            try:
                from openmm import unit
            except Exception:
                unit = None

        positions = None
        try:
            if PARMED_AVAILABLE and isinstance(ligand_struct, pmd.Structure):
                # ParmEd stores coordinates in angstroms
                coords = ligand_struct.positions
                if unit:
                    positions = [coord * unit.angstrom for coord in coords]
                else:
                    positions = coords
            elif OPENFF_AVAILABLE and hasattr(ligand_struct, 'conformers'):
                # OpenFF Molecule: conformers as numpy arrays (angstrom)
                confs = ligand_struct.conformers
                if len(confs) > 0:
                    coords = confs[0]
                    if unit:
                        positions = [tuple(coord) * unit.angstrom for coord in coords]
                    else:
                        positions = [tuple(coord) for coord in coords]
            else:
                if self.verbose:
                    print("• get_ligand_positions_openmm: unsupported ligand_struct type")
                return None

            if self.verbose:
                print(f"✓ Obtained {len(positions)} ligand positions")
            return positions
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error obtaining ligand positions: {e}")
            return None
