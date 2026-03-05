from openmm import app
import openmm.unit as unit

ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
top = app.Topology() # empty top

try:
    ff.createSystem(top, nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)
    print("SUCCESS")
except Exception as e:
    print("ERROR:", e)

