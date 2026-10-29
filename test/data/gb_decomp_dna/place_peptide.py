"""Place the peptide against the DNA strand and relax the complex in implicit solvent.

How this fixture was made (AmberTools 24.8 tleap, OpenMM 8.5):
    tleap -f leap.in            # complex/receptor/ligand .prmtop, coordinates of tleap's sequence builder
    python place_peptide.py     # writes complex_relaxed.inpcrd, stored here as complex.inpcrd
The receptor is a single DNA strand (OL15), the ligand a pentapeptide with charged termini
(ff14SB); mbondi3 radii. Used by test/unit/test_decomposition_exactness.py for
GBn2's nucleic-acid parameters and Lorentz-Berthelot (Amber) van der Waals.
"""
import numpy as np
import parmed as pmd
from openmm import app, openmm, unit

c = pmd.load_file('complex.prmtop', 'complex.inpcrd')
x = c.coordinates.copy()
dna = np.array([a.idx for a in c.atoms if a.residue.name.startswith('D')])
pep = np.array([a.idx for a in c.atoms if not a.residue.name.startswith('D')])
print('atoms', len(c.atoms), 'dna', len(dna), 'pep', len(pep))
print('dna extent', np.ptp(x[dna], axis=0), 'pep extent', np.ptp(x[pep], axis=0))
# put the peptide's centre 9 A from the DNA's centre along the axis of smallest DNA extent
ax = int(np.argmin(np.ptp(x[dna], axis=0)))
shift = x[dna].mean(0) - x[pep].mean(0)
shift[ax] += 9.0
x[pep] += shift
d = np.sqrt(((x[dna][:, None] - x[pep][None]) ** 2).sum(-1))
print('min dna-pep distance before relax', d.min())
system = c.createSystem(implicitSolvent=app.OBC2, nonbondedMethod=app.NoCutoff, constraints=None)
# pull the two together a little with a weak centroid restraint, then relax
f = openmm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-d0)^2")
f.addGlobalParameter('k', 200.0)
f.addGlobalParameter('d0', 0.75)
f.addGroup([int(i) for i in dna])
f.addGroup([int(i) for i in pep])
f.addBond([0, 1], [])
system.addForce(f)
ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName('CPU'))
ctx.setPositions(x * 0.1)
openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 500)
ctx.setParameter('k', 0.0)
openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 300)
y = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
d = np.sqrt(((y[dna][:, None] - y[pep][None]) ** 2).sum(-1))
print('min dna-pep distance after relax', d.min(), 'pairs < 4 A', (d < 4).sum())
print('complex extent', np.ptp(y, axis=0))
c.coordinates = y
c.save('complex_relaxed.inpcrd', format='rst7', overwrite=True)
