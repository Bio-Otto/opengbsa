
import pickle
import openmm
from openmm import unit
import sys
import numpy as np

def check_pickle_charge(pkl_file):
    print(f"Loading system from {pkl_file}")
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)

    # Cache format: (system, topology, positions)
    # OR dictionary
    system = None
    
    print(f"Data type: {type(data)}")
    
    if isinstance(data, dict):
        if 'system_xml' in data:
            print("Deserializing system from XML...")
            system = openmm.XmlSerializer.deserialize(data['system_xml'])
        elif 'system' in data:
            system = data['system']
        else:
             print(f"Dict keys: {data.keys()}")
             return
    elif isinstance(data, (list, tuple)):
         system = data[0]
    else:
         # Maybe the pickle IS the system?
         if isinstance(data, openmm.System):
             system = data

    if system is None:
        print("Could not extract System object from pickle!")
        return

    print(f"System particles: {system.getNumParticles()}")
    
    # Iterate forces
    nb_force = None
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            nb_force = f
            break
            
    if nb_force is None:
        print("No NonbondedForce found!")
        return

    total_charge = 0.0
    charges = []
    
    for i in range(system.getNumParticles()):
        chg, sig, eps = nb_force.getParticleParameters(i)
        q = chg.value_in_unit(unit.elementary_charge)
        charges.append(q)
        total_charge += q
        
    charges = np.array(charges)
    
    print(f"Net Charge: {total_charge:.4f} e")
    print(f"Min Charge: {np.min(charges):.4f} e")
    print(f"Max Charge: {np.max(charges):.4f} e")
    
    # Print high charges
    high_mask = np.abs(charges) > 1.0
    if np.any(high_mask):
         print("Warning: Atoms with |q| > 1.0 found!")
         print(charges[high_mask])

if __name__ == "__main__":
    check_pickle_charge(sys.argv[1])
