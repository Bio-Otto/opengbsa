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


def _lattice_shift(delta, box_vectors, inverse):
    """Lattice vector nearest to `delta` (a vector or rows of vectors), for a
    triclinic cell with box vectors as the rows of `box_vectors`.

    Rounding the reduced (fractional) coordinates gives the right image for a
    vector much shorter than the cell (a bond); for longer vectors (centres of
    molecules) the 27 neighbouring lattice vectors are searched for the
    shortest remainder, which is the true minimum image in a skewed cell such
    as a rhombic dodecahedron.
    """
    import numpy as np
    from itertools import product
    delta = np.atleast_2d(np.asarray(delta, dtype=np.float64))
    n0 = np.round(delta @ inverse)
    offsets = np.array(list(product((-1, 0, 1), repeat=3)), dtype=np.float64)
    candidates = (n0[:, None, :] + offsets[None, :, :]) @ box_vectors        # (k, 27, 3)
    best = np.argmin(np.linalg.norm(delta[:, None, :] - candidates, axis=2), axis=1)
    return candidates[np.arange(len(delta)), best]


def reimage_molecules(traj):
    """
    Re-image a trajectory across periodic boundaries, equivalent to
    `gmx trjconv -pbc mol -center` (centered on the largest molecule(s)),
    with no external tool required.

    A raw MD trajectory (e.g. GROMACS XTC) does not guarantee that molecules
    stay whole or centered -- as the simulation progresses a protein or
    ligand can drift and get wrapped to the opposite side of the periodic
    box, which does not change the physics but corrupts every energy
    calculation done directly on those raw coordinates (see
    TrajectoryProcessor.load_and_process for the numerical evidence).

    The wrapping is implemented here in numpy, bond-by-bond per frame, and
    molecules are placed in the image nearest the anchor molecules
    (cpptraj `autoimage` semantics), rather than relying on mdtraj's
    `image_molecules()` or `gmx trjconv -pbc mol`, for reasons observed on
    real production data:

    1. **Anchor heuristic.** mdtraj's default `guess_anchor_molecules()`
       heuristic is degenerate for systems with few molecules: with <= 9
       molecules the cutoff equals the size of the largest molecule, so no
       molecule can ever qualify (observed on a 3-chain antibody/PD-1
       complex: "Could not find any anchor molecules ... more than 1849
       atoms").

    2. **Frame-to-frame continuity.** mdtraj's C-level imaging assumes the
       first frame is whole and unwraps later frames relative to previous
       ones, which leaves molecules split across the box when a molecule
       jumps more than half a box between saved frames (observed at 41 ns
       of that same complex: an intra-residue HIS CB-CG bond left 10 nm
       long). Per-frame bond-by-bond wrapping has no such assumption.

    3. **Bound-partner separation.** `trjconv -pbc mol` wraps each
       molecule's center of mass into the primary cell independently; when
       a multi-chain complex straddles a box boundary this can place a
       chain on the opposite side of the box from its bound partners
       (observed at 34/41/68 ns of that same complex, where trjconv's
       output put one chain ~9.5 nm away from the other two, corrupting
       inter-chain energies). Placing every molecule in the image nearest
       the anchors avoids this -- PROVIDED the anchors are first brought into
       a single image themselves (see 4).

    4. **Anchors split across images.** When several chains qualify as
       anchors and the complex straddles a box boundary, the anchors' mean
       position must not be computed until the anchors share one image;
       otherwise it lies between images and chains are placed half a box
       apart. An earlier version did exactly that, and on antibody/PD-1
       complexes separated the antibody from PD-1 by 25-52 A in 2.5-7% of
       frames of trajectories that were intact before re-imaging.

    5. **Triclinic cells.** GROMACS rhombic dodecahedra and truncated
       octahedra are triclinic: their box vectors are not orthogonal. Earlier
       versions read only the box lengths and imaged as if the cell were
       rectangular. Triclinic cells are now imaged with the full box vectors
       (reduced coordinates through the inverse box matrix, and the nearest of
       the neighbouring images for molecule centres); rectangular cells keep
       the original arithmetic exactly.

    Steps: rebuild a standard bond graph if the topology has none (the
    temporary PDB written by TprLoader carries no CONECT records); walk
    every molecule's bond tree to make it whole; bring the anchor molecules
    into one image (largest first, each next one nearest the centre of those
    already placed); place every other molecule in the image nearest the
    anchors' centre; then shift the anchors' center of mass to the box center
    (the `-center` part of trjconv).

    Parameters
    ----------
    traj : md.Trajectory
        Trajectory with periodic unit cell information.

    Returns
    -------
    md.Trajectory
        Re-imaged copy; the input trajectory is not modified.
    """
    import numpy as np
    from collections import deque

    top = traj.topology

    # 1) Rebuild a standard bond graph if the topology has none.
    if top.n_bonds == 0:
        try:
            top.create_standard_bonds()
        except Exception as e:
            log.warning(f"Could not rebuild standard bonds ({e}); "
                        f"molecule imaging may fail.")

    molecules = list(top.find_molecules())
    if not molecules:
        raise ValueError("Trajectory topology has no molecules (no usable "
                         "bond graph); cannot re-image.")

    box = traj.unitcell_lengths
    if box is None:
        raise ValueError("This trajectory does not define a periodic unit cell")
    angles = traj.unitcell_angles
    triclinic = angles is not None and not np.allclose(angles, 90.0, atol=1e-3)
    vectors = traj.unitcell_vectors.astype(np.float64) if triclinic else None

    result = traj[:]
    xyz = result.xyz
    n_frames = xyz.shape[0]

    def frame_shifts(f):
        """(bond shift, nearest-image shift) functions of a displacement in frame f."""
        if not triclinic:
            L = box[f]
            shift = lambda delta: np.round(delta / L) * L       # rectangular: as before
            return shift, shift
        B = vectors[f]
        B_inv = np.linalg.inv(B)
        bond = lambda delta: (np.round(delta @ B_inv) @ B).astype(xyz.dtype)
        nearest = lambda delta: _lattice_shift(delta, B, B_inv)[0].astype(xyz.dtype)
        return bond, nearest

    bonds = np.asarray([[b0.index, b1.index] for b0, b1 in top.bonds],
                       dtype=np.int64)
    mol_indices = [np.asarray([a.index for a in m], dtype=np.int64)
                   for m in molecules]

    # Anchors: largest molecule plus anything at least half its size --
    # robust for multi-chain complexes whose inter-chain links (e.g.
    # engineered disulfides) are absent from the standard bond set.
    largest = max(len(m) for m in molecules)
    anchors = [m for m in molecules if len(m) >= max(1, largest // 2)]
    anchor_indices = sorted({a.index for m in anchors for a in m})
    # Anchors as index arrays, largest first: the order in which they are
    # brought into a common image in step 2b below.
    anchor_mols = sorted((np.asarray([a.index for a in m], dtype=np.int64)
                          for m in anchors), key=len, reverse=True)
    anchor_set = set(anchor_indices)

    if len(bonds):
        # adjacency for bond-tree walking
        n_atoms = xyz.shape[1]
        adj = [[] for _ in range(n_atoms)]
        for b0, b1 in bonds:
            adj[b0].append(b1)
            adj[b1].append(b0)

        # 2) Per-frame, per-molecule wrapping (`-pbc mol`): walk each
        # molecule's bond tree from an arbitrary root and shift every atom
        # into the primary image of its already-placed bonded neighbour.
        roots = [int(mi[0]) for mi in mol_indices]
        for f in range(n_frames):
            bond_shift, image_shift = frame_shifts(f)
            for root, mi in zip(roots, mol_indices):
                visited = {root}
                queue = deque([root])
                while queue:
                    i = queue.popleft()
                    xi = xyz[f, i]
                    for j in adj[i]:
                        if j not in visited:
                            visited.add(j)
                            delta = xyz[f, j] - xi
                            xyz[f, j] -= bond_shift(delta)
                            queue.append(j)

            # 2b) Bring the ANCHORS into one image first. When a multi-chain
            # complex straddles a box boundary, its chains sit in different
            # periodic images; averaging them as they are gives a centre that
            # belongs to none of them, and placing each chain nearest that
            # point can put a bound chain a box length from its partners.
            # This needs a complex that is large relative to the box (chain
            # centres spread over ~L/2, e.g. an antibody Fv/PD-1 complex in a
            # cubic box with 1 nm padding); compact complexes are unaffected,
            # which is why it went unnoticed. Regression test:
            # test/unit/test_reimage_multichain.py.
            # This happened: on antibody/PD-1 complexes where all three chains
            # qualify as anchors, the previous version (which averaged the
            # anchors before placing them) separated the antibody from PD-1 by
            # 25-52 A in 20/800 and 55/800 frames of trajectories that were
            # intact in the raw coordinates -- zero interface energy, no
            # warning; found by cross-validation against Amber MMPBSA.py.
            # So: keep the largest anchor where it is, then place each further
            # anchor in the image nearest the centre of those already placed.
            placed = [anchor_mols[0]]
            for mi in anchor_mols[1:]:
                ref = xyz[f, np.concatenate(placed)].mean(axis=0)
                com = xyz[f, mi].mean(axis=0)
                xyz[f, mi] -= image_shift(com - ref)
                placed.append(mi)

            # 2c) Place every remaining molecule in the image NEAREST the
            # (now coherent) anchor centre of mass (cpptraj `autoimage`
            # behaviour). Wrapping each molecule's COM into the primary cell
            # instead (what `trjconv -pbc mol` does) can separate bound
            # partners: when a multi-chain complex straddles a box boundary,
            # chains are wrapped into opposite sides of the box and
            # inter-chain distances become box-length large (observed at
            # 34/41/68 ns of a 3-chain antibody/PD-1 trajectory, where
            # trjconv's output put one chain ~9.5 nm from the other two).
            anchor_com = xyz[f, anchor_indices].mean(axis=0)
            for mi in mol_indices:
                if int(mi[0]) in anchor_set:
                    continue
                com = xyz[f, mi].mean(axis=0)
                xyz[f, mi] -= image_shift(com - anchor_com)

    # 3) Center anchor COM in the box (the `-center` part of trjconv); the
    # centre of a triclinic cell is half the sum of its box vectors.
    com = xyz[:, anchor_indices, :].mean(axis=1)
    centre = 0.5 * vectors.sum(axis=1) if triclinic else box / 2.0
    xyz += (centre - com).astype(xyz.dtype)[:, np.newaxis, :]

    return result


class TrajectoryProcessor:
    """
    Handles loading, slicing, and processing (stripping) of trajectories.
    """

    @staticmethod
    def resolve_mdtraj_loadable_path(trajectory_file, alias_dir=None):
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

        The alias is created in `alias_dir` (the run folder), or in a new
        temporary folder when none is given -- never next to the input.
        """
        if not str(trajectory_file).lower().endswith('.trj'):
            return str(trajectory_file)
        if alias_dir is None:
            import tempfile
            alias_dir = tempfile.mkdtemp(prefix='opengbsa_trj_')
        alias_dir = Path(alias_dir)
        alias_dir.mkdir(parents=True, exist_ok=True)
        alias_path = alias_dir / (Path(str(trajectory_file)).name[:-4] + '.mdcrd')
        target = Path(trajectory_file).resolve()
        if alias_path.is_symlink() and alias_path.resolve() != target:
            alias_path.unlink()
        if not alias_path.exists():
            try:
                alias_path.symlink_to(target)
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
                traj = reimage_molecules(traj)
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
