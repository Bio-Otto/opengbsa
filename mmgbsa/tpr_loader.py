"""
Native TPR Parser logic bypassing ParmEd's outdated parser.
Works for GROMACS 2023+ (Version 137+).
Uses the TprParser C++ Python bindings wrapper.
"""
import parmed as pmd
import numpy as np

import mdtraj as md

def load_tpr_as_parmed(tpr_path: str, xtc_path: str = None) -> pmd.Structure:
    from TprParser.TprReader import TprReader
    reader = TprReader(tpr_path)
    
    struct = pmd.Structure()
    
    atomnames = reader.get_name('atom')
    resnames = reader.get_name('res')
    resids = reader.get_ivector('resid')
    masses = reader.get_mq('m')
    charges = reader.get_mq('q')
    coords = reader.get_xvf('x') # nm
    n_atoms = len(atomnames)

    # --- Virtual site and Dummy Atom Detection ---
    vsites_raw = None
    try:
        vsites_raw = reader.get_vsites()
    except Exception:
        pass

    vsite_atom_indices_1based = set()
    if vsites_raw is not None and len(vsites_raw) > 0:
        for vs in vsites_raw:
            vsite_atom_indices_1based.add(int(vs[1]))

    mass_zero_indices = set(np.where(np.array(masses, dtype=float) == 0.0)[0] + 1)
    all_dummy_1based = vsite_atom_indices_1based | mass_zero_indices
    dummy_indices = np.array(sorted(all_dummy_1based), dtype=int) - 1
    real_indices = np.array([i for i in range(n_atoms) if i not in set(dummy_indices)], dtype=int)
    
    use_indices = np.arange(n_atoms)
    
    # --- XTC vs TPR Atom Count Validation (Dynamically Exclude Dummies) ---
    if xtc_path is not None:
        try:
            import mdtraj as md
            import openmm.app as app
            
            top = app.Topology()
            chain = top.addChain()
            prev_resid = None
            residue = None
            for i in real_indices:
                rid = int(resids[i])
                if rid != prev_resid:
                    residue = top.addResidue(resnames[i], chain)
                    prev_resid = rid
                m = float(masses[i])
                element = app.Element.getByMass(m) if m > 0 else None
                top.addAtom(atomnames[i], element, residue)
                
            mdtraj_top = md.Topology.from_openmm(top)
            first_chunk = next(md.iterload(xtc_path, top=mdtraj_top, chunk=1))
            n_atoms_xtc = first_chunk.n_atoms
            
            if n_atoms_xtc == len(real_indices):
                use_indices = real_indices
                print(f"[TprLoader] XTC has {n_atoms_xtc} atoms. Excluding {len(dummy_indices)} dummy atoms from topology.")
            elif n_atoms_xtc == n_atoms:
                use_indices = np.arange(n_atoms)
                print(f"[TprLoader] XTC has {n_atoms_xtc} atoms. Including all dummy atoms.")
            else:
                print(f"[TprLoader] WARNING: XTC atom count ({n_atoms_xtc}) matches neither full ({n_atoms}) nor real ({len(real_indices)}) counts. Defaulting to all.")
        except Exception as e:
            print(f"[TprLoader] WARNING: Failed to verify XTC atom count: {e}. Defaulting to all atoms.")
            use_indices = np.arange(n_atoms)
    
    use_set = set(use_indices)
    
    # Map from TPR index to ParmEd index
    tpr_to_pmd = {}
    
    for pmd_idx, i in enumerate(use_indices):
        aname = atomnames[i]
        rname = resnames[i]
        rid = resids[i]
        charge = charges[i]
        mass = masses[i]
        
        # Approximate atomic number from mass since GROMACS abstracts it
        atomic_number = 0
        
        if 1.0 <= mass < 3.0: atomic_number = 1   # H
        elif 11.0 <= mass < 14.0: atomic_number = 6 # C
        elif 14.0 <= mass < 15.0: atomic_number = 7 # N
        elif 15.0 <= mass < 18.0: atomic_number = 8 # O
        elif 30.0 <= mass < 32.0: atomic_number = 15 # P
        elif 32.0 <= mass < 35.0: atomic_number = 16 # S
        elif 35.0 <= mass < 37.0: atomic_number = 17 # Cl
        elif 18.0 <= mass < 20.0: atomic_number = 9 # F
        elif 22.0 <= mass < 24.0: atomic_number = 11 # Na
        elif 24.0 <= mass < 26.0: atomic_number = 12 # Mg
        elif 39.0 <= mass < 41.0: atomic_number = 19 # K
        elif 40.0 <= mass < 42.0: atomic_number = 20 # Ca
        else:
            # Fallback to name heuristic
            elem = aname[0:1] if aname[0:1] in "C N O S P F I H K" else aname[0:2].upper()
            try:
                import parmed.periodic_table as pt
                atomic_number = pt.AtomicNum.get(elem, 0)
            except:
                pass
        
        atom = pmd.Atom(name=aname, type=aname, charge=charge, mass=mass, atomic_number=atomic_number)
        struct.add_atom(atom, resname=rname, resnum=rid)
        tpr_to_pmd[i] = pmd_idx
        
    bonds = reader.get_bonded('bonds')
    dummy_bond = pmd.BondType(100.0, 1.0)
    struct.bond_types.append(dummy_bond)
    dummy_bond.list = struct.bond_types
    
    for b in bonds:
        idx1 = int(b[0]) - 1
        idx2 = int(b[1]) - 1
        if idx1 in tpr_to_pmd and idx2 in tpr_to_pmd:
            a1 = struct.atoms[tpr_to_pmd[idx1]]
            a2 = struct.atoms[tpr_to_pmd[idx2]]
            bond = pmd.Bond(a1, a2)
            bond.type = dummy_bond
            struct.bonds.append(bond)
        
    # Nonbonded LJ parameters
    lj = reader.get_nonbonded('lj')
    atom_types = {}
    for pmd_idx, i in enumerate(use_indices):
        atom = struct.atoms[pmd_idx]
        lj_data = lj[i]
        type_idx = int(lj_data[0])
        sigma_nm = lj_data[1]
        epsilon_kj = lj_data[2]
        
        sigma_ang = sigma_nm * 10.0
        epsilon_kcal = epsilon_kj / 4.184
        
        if type_idx not in atom_types:
            atype = pmd.AtomType(name=f"t{type_idx}", number=type_idx, mass=atom.mass)
            atype.sigma = sigma_ang
            atype.epsilon = epsilon_kcal
            atom_types[type_idx] = atype
            
        atom.atom_type = atom_types[type_idx]
        atom.sigma = sigma_ang
        atom.epsilon = epsilon_kcal

    subset_coords = [coords[i] for i in use_indices]
    coords_ang = np.array(subset_coords) * 10.0
    struct.coordinates = coords_ang
    
    # Needs a box so trajectory files don't fail length checks if they have unitcells
    try:
        box = reader.get_xvf('box')
        # TprParser box is a 3x3 matrix in nm
        if len(box) == 3:
            # simple assumed orthogonal
            v1, v2, v3 = box[0][0], box[1][1], box[2][2]
            struct.box = [v1*10.0, v2*10.0, v3*10.0, 90.0, 90.0, 90.0]
    except:
        pass

    return struct
