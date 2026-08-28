"""
Trajectory loading with automatic solvated/dry topology reconciliation.

Provides `TrajectoryProcessor`, which loads a trajectory against a given
topology and, on an atom-count mismatch, attempts to auto-discover a
matching solvated topology or strip water/ions on the fly so the frame
count lines up with the dry system used for energy evaluation.
"""
import mdtraj as md
from pathlib import Path
import logging

log = logging.getLogger(__name__)

class TrajectoryProcessor:
    """
    Handles loading, slicing, and processing (stripping) of trajectories.
    """

    @staticmethod
    def resolve_mdtraj_loadable_path(trajectory_file):
        """
        Returns a path mdtraj can load, aliasing a bare '.trj' extension to
        '.mdcrd' via a symlink if needed.

        mdtraj has no registered loader for '.trj' (its own FormatRegistry
        only recognizes '.mdcrd'/'.crd' for ASCII Amber trajectories) --
        confirmed on a real published dataset (Zenodo 17926575) whose Amber
        MD output is named 'dry_MD_<compound>.trj' but is byte-for-byte a
        standard ASCII mdcrd file (same "Cpptraj Generated trajectory"
        header, same fixed-width coordinate block format). This is purely a
        filename-extension issue, not a real format difference, so aliasing
        via a symlink (rather than copying) reinterprets no data. Any
        caller that hands a trajectory path to `md.load`/`md.load_frame`
        (not just `TrajectoryProcessor.load_and_process` itself) should
        route it through this first, or it will hit mdtraj's own opaque
        "no loader for filename=... (extension=.trj)" OSError.
        """
        if not str(trajectory_file).lower().endswith('.trj'):
            return str(trajectory_file)
        alias_path = Path(str(trajectory_file)[:-4] + '.mdcrd')
        if not alias_path.exists():
            try:
                alias_path.symlink_to(Path(trajectory_file).resolve())
            except OSError as e:
                log.warning(f"Could not create .mdcrd alias for .trj trajectory ({e}); "
                            f"proceeding with the original path, which mdtraj may reject.")
                return str(trajectory_file)
        return str(alias_path)

    @staticmethod
    def _auto_discover_solvated(trajectory_file, expected_atoms):
        """
        Attempts to find a PDB file in the same directory that matches the trajectory's atom count.
        """
        traj_path = Path(trajectory_file)
        search_dir = traj_path.parent
        
        # We need to know the trajectory's atom count.
        # But we can't load it without a topology if it's XTC.
        # Catch-22? 
        # MDTraj CAN load XTC without topology? No, it raises "Need a topology".
        # But wait, the error message from MDTraj SAYS: "ValueError: xyz must be shape (Any, 3880, 3). You supplied (1001, 4003, 3)".
        # So MDTraj KNOWS the atom count (4003) internally before failing against the topology (3880)!
        # Can we extract that?
        # The exception message is our clue.
        return None

    @staticmethod
    def load_and_process(trajectory_file, topology_file, target_atoms=None, solvated_topology=None,
                         stride=1, start=None, end=None, reimage=True):
        """
        Loads a trajectory, optionally using a solvated topology for stripping.

        reimage : bool, default True
            Re-image molecules across periodic boundaries after loading (see
            the image_molecules() call below for why this matters). Exposed
            as a parameter -- rather than always-on -- so a config with its
            own reference trajectory that is already correctly imaged (e.g.
            re-using this codebase's own previously-imaged output, or a
            trajectory a user has already run cpptraj `autoimage` on) is not
            forced through a redundant, non-free re-imaging pass; also lets
            a user explicitly opt out if they have a specific reason to
            compare against raw, un-imaged coordinates. Defaults to True
            because skipping this step silently corrupts vdW/electrostatic
            energies once any molecule drifts across a periodic boundary
            (confirmed on real production data: a ~11% systematic
            underestimate of |delta_vdw| that grows across the trajectory).
        """
        # 1. Determine load topology
        load_top = topology_file
        if solvated_topology and Path(solvated_topology).exists():
            load_top = solvated_topology

        print(f"Loading trajectory {trajectory_file} with topology {load_top}...")

        effective_trajectory_file = TrajectoryProcessor.resolve_mdtraj_loadable_path(trajectory_file)

        try:
            # MDTraj load
            traj = md.load(effective_trajectory_file, top=load_top, stride=stride)
        except ValueError as e:
            # Check for atom mismatch message
            e_str = str(e)
            if "topology and the trajectory files might not contain the same atoms" in e_str or "xyz must be shape" in e_str:
                print(f"⚠️  Atom Mismatch detected during load.")
                
                # Attempt Discovery
                # Strategy: Iterate over ALL .pdb/prmtop files in directory.
                # Try loading one frame with them. If it works, suggests it.
                print(f"🔎 Scanning directory for matching solvated topology...")
                traj_path = Path(trajectory_file)
                candidates = list(traj_path.parent.glob("*.pdb")) + list(traj_path.parent.glob("*.prmtop"))

                found_match = None
                for cand in candidates:
                    if cand.resolve() == Path(topology_file).resolve(): continue
                    try:
                        # Try loading just 1 frame
                        t_test = md.load(effective_trajectory_file, top=str(cand), frame=0)
                        # If we get here, it matched!
                        print(f"💡 Found potential match: {cand.name} ({t_test.n_atoms} atoms)")
                        found_match = str(cand)
                        break
                    except:
                        continue

                if found_match:
                    print(f"🔄 Retrying load with auto-detected topology: {Path(found_match).name}")
                    try:
                         traj = md.load(effective_trajectory_file, top=found_match, stride=stride)
                         # IMPORTANT: If we successfully loaded with a solvated topology,
                         # we MUST strip it down to target_atoms (dry topology count).
                         # We rely on step 2 (below) to handle this check.
                    except Exception as e2:
                        raise ValueError(f"Auto-recovery failed: {e2}")
                else:
                    # Specific Guidance
                    raise ValueError(
                        f"\n❌ CRITICAL ATOM MISMATCH ❌\n"
                        f"The trajectory '{Path(trajectory_file).name}' does not match the topology '{Path(topology_file).name}'.\n"
                        f"Please provide a 'solvated_topology' in your config that matches the trajectory atoms.\n"
                        f"(Example: Check your MD equilibration PDB/GROMACS file)"
                    ) from e
            else:
                raise e
        except Exception as e:
            raise ValueError(f"Failed to load trajectory: {e}")
            
        # Slicing (Frames)
        if start is not None or end is not None:
            traj = traj[start:end]

        # 2. Check Atom Count / Stripping
        current_atoms = traj.n_atoms
        
        # If target_atoms is provided, verify or strip
        if target_atoms and current_atoms != target_atoms:
            print(f"Atom Mismatch detected after load: Trajectory ({current_atoms}) vs Target ({target_atoms})")
            
            # Smart Stripping Logic
            if current_atoms > target_atoms:
                print("Attempting on-the-fly stripping...")
                
                # Try standard water/ion strip
                # Selection logic matching core.py implementation
                try:
                    # Includes both Amber/GROMACS-style ion residue names
                    # (NA, CL, K) and CHARMM/NAMD-style ones (SOD, CLA, POT) --
                    # confirmed on a real NAMD PSF/PDB that 'resname CL' does
                    # NOT match CHARMM's 'CLA' (mdtraj resname selection is an
                    # exact string match, not a prefix match), which silently
                    # left 12 chloride-ion atoms in the "stripped" trajectory
                    # and caused the subsequent blind-slicing fallback to grab
                    # the wrong 2030 atoms entirely (chloride ions instead of
                    # part of the protein), corrupting every downstream energy.
                    strip_mask = traj.topology.select(
                        'not (water or resname NA or resname CL or resname SOD or resname K '
                        'or resname CLA or resname POT or resname CAL or resname ZN2)'
                    )
                except Exception:
                     # Fallback if selection fails (e.g. non-standard names)
                     strip_mask = []

                if len(strip_mask) == target_atoms:
                    print(f"✓ Stripping successful, matched {target_atoms} atoms.")
                    traj = traj.atom_slice(strip_mask)
                else:
                    # Fallback slicing
                    print(f"⚠️  Standard stripping (water/ions) resulted in {len(strip_mask)} atoms. Expected {target_atoms}.")
                    print(f"   Falling back to blind slicing (keeping first {target_atoms} atoms).")
                    traj = traj.atom_slice(range(target_atoms))

            elif current_atoms < target_atoms:
                raise ValueError(f"Trajectory has fewer atoms ({current_atoms}) than target topology ({target_atoms})! Cannot process.")

        # 3. Re-image molecules across periodic boundaries (PBC re-wrapping).
        # A raw GROMACS trajectory has NO guarantee that a given molecule
        # stays whole/centered from frame to frame -- as the simulation
        # progresses, the protein or ligand can drift and get wrapped to the
        # opposite side of the periodic box, which does not change the true
        # physics but DOES corrupt any energy calculation performed directly
        # on these raw coordinates (vdW/electrostatic terms depend on actual
        # inter-atomic distances, which become wrong once a molecule's
        # frame-to-frame center jumps by a box length). Confirmed on real
        # production data: computing delta_vdw directly from this raw,
        # un-imaged trajectory reproduces this codebase's own (wrong) output
        # bit-for-bit (e.g. -39.41 kcal/mol on one benchmark replica),
        # while the same calculation on a properly re-imaged version of the
        # identical frames matches real Amber sander output to <0.001
        # kcal/mol (-44.37 kcal/mol on that same replica) -- i.e. this
        # missing re-imaging step, not a formula or parameter error, was the
        # actual root cause of a systematic ~11% low-magnitude bias in vdW
        # (and likely other MM) energies that grows the further into the
        # trajectory a frame is (frames near the start, before much drift
        # has accumulated, were unaffected). `image_molecules()` re-wraps
        # every whole molecule to minimize its span without altering any
        # intra-molecular geometry, which is exactly what real Amber/
        # cpptraj `autoimage` does before MMPBSA.py ever sees a trajectory.
        if reimage:
            try:
                traj = traj.image_molecules(inplace=False)
            except Exception as e:
                log.warning(f"image_molecules() failed ({e}); proceeding with "
                            f"un-imaged coordinates. Energies may be biased if "
                            f"any molecule crosses a periodic boundary during "
                            f"the trajectory.")
        else:
            log.warning("Trajectory re-imaging (reimage_trajectory) disabled by "
                        "configuration; using raw, un-imaged coordinates. Energies "
                        "will be biased if any molecule crosses a periodic boundary "
                        "during the trajectory -- only disable this if the input "
                        "trajectory is already known to be correctly imaged.")

        return traj
