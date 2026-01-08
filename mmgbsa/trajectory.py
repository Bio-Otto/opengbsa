import mdtraj as md
from pathlib import Path
import logging

log = logging.getLogger(__name__)

class TrajectoryProcessor:
    """
    Handles loading, slicing, and processing (stripping) of trajectories.
    """

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
                         stride=1, start=None, end=None):
        """
        Loads a trajectory, optionally using a solvated topology for stripping.
        """
        # 1. Determine load topology
        load_top = topology_file
        if solvated_topology and Path(solvated_topology).exists():
            load_top = solvated_topology
            
        print(f"Loading trajectory {trajectory_file} with topology {load_top}...")
            
        try:
            # MDTraj load
            traj = md.load(trajectory_file, top=load_top, stride=stride)
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
                        t_test = md.load(trajectory_file, top=str(cand), frame=0)
                        # If we get here, it matched!
                        print(f"💡 Found potential match: {cand.name} ({t_test.n_atoms} atoms)")
                        found_match = str(cand)
                        break
                    except:
                        continue
                
                if found_match:
                    print(f"🔄 Retrying load with auto-detected topology: {Path(found_match).name}")
                    try:
                         traj = md.load(trajectory_file, top=found_match, stride=stride)
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
                    strip_mask = traj.topology.select('not (water or resname NA or resname CL or resname SOD or resname K)')
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
                
        return traj
