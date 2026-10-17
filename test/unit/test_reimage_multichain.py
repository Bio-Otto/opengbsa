"""Regression test: reimage_molecules() must keep a multi-chain complex together
when the complex straddles a periodic box boundary.

Background (2026-10-11): for an antibody Fv/PD-1 complex all three chains are
"anchors" (each is at least half the size of the largest). The function used to
average the anchors' coordinates BEFORE bringing them into a common periodic
image; when the complex sat across a box face, that average landed mid-box and
the chains were placed half a box apart -- zero interface energy in 2.5-7% of
frames, with no warning. Found by cross-validation against Amber MMPBSA.py.

The condition that triggers it is a complex that is LARGE relative to the box:
the PD-1 systems used a cubic box with 1 nm padding (L ~ 10 nm), so the bound
chains' centres spread over ~5 nm. Worked example, chain centres at x = -2.0,
0.2, 3.5 nm (true, unwrapped): the first wraps to 8.0; the anchor mean is 3.9;
8.0 is only 4.1 nm from it, so the old code left that chain at 8.0 -- 7.8 nm
from its partner instead of 2.2. A compact cluster (chains within ~1 nm) does
NOT trigger it, so the fixture below deliberately uses the realistic spread.
"""
import mdtraj as md
import numpy as np

from mmgbsa.trajectory import reimage_molecules

L = 10.0                      # cubic box edge, nm (PD-1 systems: 10.06-10.91 nm)
TRUE_X = (-2.0, 0.2, 3.5)     # chain centres in the bound complex, unwrapped
N_PER_CHAIN = 12
SPACING = 0.15                # nm between bonded atoms


def _three_chain_complex(n_frames=3):
    top = md.Topology()
    for _ in range(3):
        chain = top.add_chain()
        prev = None
        for _ in range(N_PER_CHAIN):
            res = top.add_residue("ALA", chain)
            atom = top.add_atom("CA", md.element.carbon, res)
            if prev is not None:
                top.add_bond(prev, atom)
            prev = atom

    # Bound geometry: chains side by side along x (contacts 0-1 and 1-2), each a
    # short rod along y; then wrap into the primary cell as MD output does, which
    # moves chain 0 (true x = -2.0) to x = 8.0, across the box from its partner.
    ys = np.arange(N_PER_CHAIN) * SPACING + 1.0
    frame = []
    for x in TRUE_X:
        for y in ys:
            frame.append([x, y, 2.0])
    xyz = np.array(frame, dtype=np.float32)
    xyz[:, 0] %= L
    xyz = np.repeat(xyz[None], n_frames, axis=0)
    cell = np.repeat(np.array([[L, L, L]], dtype=np.float32), n_frames, axis=0)
    angles = np.repeat(np.array([[90.0, 90.0, 90.0]], dtype=np.float32), n_frames, axis=0)
    return md.Trajectory(xyz, top, unitcell_lengths=cell, unitcell_angles=angles)


def _chain_coms(traj):
    out = []
    for chain in traj.topology.chains:
        idx = [a.index for a in chain.atoms]
        out.append(traj.xyz[:, idx, :].mean(axis=1))
    return np.stack(out, axis=1)            # (frames, chains, 3)


def test_raw_coordinates_are_split_across_the_box():
    """Sanity check on the fixture: in wrapped coordinates chain 0 sits ~L-2.2 nm
    from its bound partner, chain 1."""
    com = _chain_coms(_three_chain_complex())
    assert np.abs(com[:, 1, 0] - com[:, 0, 0]).max() > L / 2


def test_reimage_restores_the_bound_geometry():
    """All three chains are anchors (equal size) -- the case that used to fail.
    After imaging, the chains' relative positions must be the bound ones
    (0-1: 2.2 nm, 1-2: 3.3 nm, 0-2: 5.5 nm along x), up to a common shift."""
    imaged = reimage_molecules(_three_chain_complex())
    com = _chain_coms(imaged)[..., 0]                  # x components, (frames, chains)
    rel = com - com[:, [0]]
    expected = np.array(TRUE_X) - TRUE_X[0]
    np.testing.assert_allclose(rel, np.broadcast_to(expected, rel.shape), atol=1e-4,
                               err_msg="bound chains were placed in different periodic images")


def test_reimage_does_not_modify_input():
    traj = _three_chain_complex()
    before = traj.xyz.copy()
    reimage_molecules(traj)
    np.testing.assert_array_equal(traj.xyz, before)
