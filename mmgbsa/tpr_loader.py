"""
Native TPR Parser logic bypassing ParmEd's outdated parser.
Works for GROMACS 2023+ (Version 137+).
Uses the TprParser C++ Python bindings wrapper.
"""
import parmed as pmd
import numpy as np

import mdtraj as md

def load_tpr_as_parmed(tpr_path: str, xtc_path: str = None) -> pmd.Structure:
    """
    Build a ParmEd Structure directly from a GROMACS .tpr file's binary run-input
    data (atoms, charges, masses, LJ parameters, bonded force field terms, box),
    bypassing ParmEd's own (older) TPR reader.

    Parameters
    ----------
    tpr_path : str
        Path to the GROMACS .tpr file.
    xtc_path : str, optional
        Path to a matching trajectory. If given, used only to disambiguate
        whether virtual-site/dummy atoms present in the TPR should be included
        in the returned Structure, by comparing atom counts against the first
        trajectory frame.

    Returns
    -------
    parmed.Structure
        Atoms, LJ nonbonded parameters, harmonic bonds, harmonic/Urey-Bradley
        angles, proper dihedrals (GROMACS funct 1/9), and impropers (funct 2
        harmonic / funct 4 periodic) are parsed and added with real force
        field parameters extracted from the TPR.

    Known limitations
    -----------------
    - **CMAP (CHARMM backbone grid correction) is not parsed.** If the input
      topology uses a CHARMM-family force field with CMAP terms (as all KCX/
      OXA-type systems in this package do), those 2D phi/psi grid correction
      energies are silently absent from the returned Structure and from any
      OpenMM System built from it. A warning is printed when CMAP terms are
      detected in the TPR but this is not a substitute for actually applying
      them. Prefer `.prmtop` Native mode (which goes through ParmEd's AMBER
      parser, unaffected by this limitation) when CMAP-accurate absolute
      internal energies matter; for ligand binding free energy differences
      (complex - receptor - ligand), CMAP contributions to the unchanged
      protein backbone largely cancel, but this has not been verified for the
      per-residue decomposition of frames with non-negligible backbone motion.
    - Ryckaert-Bellemans dihedrals (GROMACS funct 3) are not parsed (CHARMM
      force fields, the primary target of this package, do not produce them).
    - Raw atomic coordinates are read as stored in the TPR without periodic-
      boundary molecular re-wrapping: any solvent/ion molecule split across a
      box edge appears with a stretched bond in this Structure's own
      coordinates. In the standard pipeline this is inconsequential because
      solvent is stripped from the Structure before any OpenMM System/energy
      is built from it, and per-frame trajectory coordinates for the actual
      GBSA calculation come from the (properly-wrapped) .xtc file, not from
      this function's embedded coordinates -- but any other use of this
      function's `.coordinates`/`.positions` directly (e.g. treating the
      returned Structure as a general-purpose PDB-equivalent single frame)
      should not assume solvent geometry is whole/unwrapped-safe.
    - The GROMACS box is read as a general triclinic cell (lengths + angles
      derived from the full 3x3 vector matrix), but any box vector convention
      other than GROMACS' row-major lower-triangular storage is not handled.
    """
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
        
    # --- Bonded force field parameters (bonds/angles/dihedrals/impropers) ---
    # GROMACS reports A-state (and, for FEP topologies, B-state) parameters per
    # entry, tagged with a GROMACS "function type" integer. Unit conversions:
    #   length: nm -> Angstrom (x10)
    #   bond/UB force constant: kJ/mol/nm^2 -> kcal/mol/Angstrom^2 (/4.184/100)
    #   angle force constant: kJ/mol/rad^2 -> kcal/mol/rad^2 (/4.184)
    #   dihedral/improper force constant: kJ/mol -> kcal/mol (/4.184)
    KJ_NM2_TO_KCAL_A2 = 1.0 / (4.184 * 100.0)
    KJ_TO_KCAL = 1.0 / 4.184
    bond_type_cache = {}

    def _get_bond_type(b0_nm, kb_kjmolnm2):
        """Harmonic bond/Urey-Bradley type, cached by (b0, kb) to avoid duplicate BondType objects."""
        req = b0_nm * 10.0
        k = kb_kjmolnm2 * KJ_NM2_TO_KCAL_A2
        key = (round(req, 6), round(k, 6))
        bt = bond_type_cache.get(key)
        if bt is None:
            bt = pmd.BondType(k, req)
            struct.bond_types.append(bt)
            bt.list = struct.bond_types
            bond_type_cache[key] = bt
        return bt

    bonds = reader.get_bonded('bonds')
    n_bonds_no_k = 0
    for b in bonds:
        idx1 = int(b[0]) - 1
        idx2 = int(b[1]) - 1
        if idx1 not in tpr_to_pmd or idx2 not in tpr_to_pmd:
            continue
        functype = int(b[2])
        params = b[3:]
        if functype != 1:
            # Only the standard harmonic bond potential (funct 1) is handled;
            # anything else falls back to the funct-1 field layout best-effort.
            pass
        if len(params) >= 4:
            b0_nm, kb = float(params[0]), float(params[1])
        elif len(params) >= 2:
            # GROMACS omits the force constant for bonds whose length is fully
            # constrained (e.g. X-H bonds under LINCS/SETTLE at MD time). Since
            # this package evaluates single-point energies with OpenMM (no
            # constraint solver in play for energy evaluation), such bonds are
            # given a very stiff, physically representative constant (the
            # typical AMBER/CHARMM X-H stretch force constant, ~340 kcal/mol/A^2)
            # rather than an arbitrary placeholder, since they still contribute
            # a (small, near-equilibrium) bonded energy term.
            b0_nm, kb = float(params[0]), 340.0 / KJ_NM2_TO_KCAL_A2
            n_bonds_no_k += 1
        else:
            continue
        a1 = struct.atoms[tpr_to_pmd[idx1]]
        a2 = struct.atoms[tpr_to_pmd[idx2]]
        bond = pmd.Bond(a1, a2)
        bond.type = _get_bond_type(b0_nm, kb)
        struct.bonds.append(bond)
    if n_bonds_no_k:
        print(f"[TprLoader] {n_bonds_no_k} constrained bond(s) had no GROMACS force constant "
              f"(LINCS/SETTLE-constrained); assigned a representative stiff harmonic constant "
              f"for single-point energy evaluation.")

    # --- Angles (harmonic, funct 1; Urey-Bradley, funct 5) ---
    angle_type_cache = {}
    ub_type_cache = {}

    def _get_angle_type(theta0_deg, ktheta_kjmolrad2):
        k = ktheta_kjmolrad2 * KJ_TO_KCAL
        key = (round(theta0_deg, 6), round(k, 6))
        at = angle_type_cache.get(key)
        if at is None:
            at = pmd.AngleType(k, theta0_deg)
            struct.angle_types.append(at)
            at.list = struct.angle_types
            angle_type_cache[key] = at
        return at

    angles = reader.get_bonded('angles')
    for ang in angles:
        i1, i2, i3 = int(ang[0]) - 1, int(ang[1]) - 1, int(ang[2]) - 1
        if i1 not in tpr_to_pmd or i2 not in tpr_to_pmd or i3 not in tpr_to_pmd:
            continue
        functype = int(ang[3])
        params = ang[4:]
        if len(params) < 2:
            continue
        theta0, ktheta = float(params[0]), float(params[1])
        a1 = struct.atoms[tpr_to_pmd[i1]]
        a2 = struct.atoms[tpr_to_pmd[i2]]
        a3 = struct.atoms[tpr_to_pmd[i3]]
        angle = pmd.Angle(a1, a2, a3)
        angle.type = _get_angle_type(theta0, ktheta)
        struct.angles.append(angle)

        if functype == 5 and len(params) >= 4:
            # Urey-Bradley 1-3 term (CHARMM): r13 (nm), kUB (kJ/mol/nm^2)
            r13_nm, kub = float(params[2]), float(params[3])
            if kub != 0.0:
                key = (round(r13_nm * 10.0, 6), round(kub * KJ_NM2_TO_KCAL_A2, 6))
                ubt = ub_type_cache.get(key)
                if ubt is None:
                    ubt = _get_bond_type(r13_nm, kub)
                    ub_type_cache[key] = ubt
                ub = pmd.UreyBradley(a1, a3)
                ub.type = ubt
                struct.urey_bradleys.append(ub)

    # --- Proper dihedrals (funct 1 = single term, funct 9 = multi-term periodic) ---
    dihedral_type_cache = {}

    def _get_dihedral_type(phi0_deg, kphi_kjmol, mult):
        k = kphi_kjmol * KJ_TO_KCAL
        key = (round(phi0_deg, 6), round(k, 6), int(mult))
        dt = dihedral_type_cache.get(key)
        if dt is None:
            dt = pmd.DihedralType(k, int(mult) if mult else 1, phi0_deg)
            struct.dihedral_types.append(dt)
            dt.list = struct.dihedral_types
            dihedral_type_cache[key] = dt
        return dt

    dihedrals = reader.get_bonded('dihedrals')
    for dih in dihedrals:
        i1, i2, i3, i4 = (int(dih[0]) - 1, int(dih[1]) - 1, int(dih[2]) - 1, int(dih[3]) - 1)
        if any(ix not in tpr_to_pmd for ix in (i1, i2, i3, i4)):
            continue
        functype = int(dih[4])
        params = dih[5:]
        if functype in (1, 9) and len(params) >= 3:
            phi0, kphi, mult = float(params[0]), float(params[1]), params[2]
            a1 = struct.atoms[tpr_to_pmd[i1]]
            a2 = struct.atoms[tpr_to_pmd[i2]]
            a3 = struct.atoms[tpr_to_pmd[i3]]
            a4 = struct.atoms[tpr_to_pmd[i4]]
            d = pmd.Dihedral(a1, a2, a3, a4, improper=False)
            d.type = _get_dihedral_type(phi0, kphi, mult)
            struct.dihedrals.append(d)
        # funct 3 (Ryckaert-Bellemans) is not expected from CHARMM-derived
        # inputs but is not handled here; such dihedrals are silently skipped
        # rather than mis-parameterized.

    # --- Impropers (funct 2 = harmonic improper; funct 4 = periodic improper) ---
    improper_type_cache = {}

    def _get_improper_type(psi0_deg, kpsi_kjmolrad2):
        k = kpsi_kjmolrad2 * KJ_TO_KCAL
        key = (round(psi0_deg, 6), round(k, 6))
        it = improper_type_cache.get(key)
        if it is None:
            it = pmd.ImproperType(k, psi0_deg)
            struct.improper_types.append(it)
            it.list = struct.improper_types
            improper_type_cache[key] = it
        return it

    try:
        impropers = reader.get_bonded('impropers')
    except Exception:
        impropers = []
    for imp in impropers:
        i1, i2, i3, i4 = (int(imp[0]) - 1, int(imp[1]) - 1, int(imp[2]) - 1, int(imp[3]) - 1)
        if any(ix not in tpr_to_pmd for ix in (i1, i2, i3, i4)):
            continue
        functype = int(imp[4])
        params = imp[5:]
        a1 = struct.atoms[tpr_to_pmd[i1]]
        a2 = struct.atoms[tpr_to_pmd[i2]]
        a3 = struct.atoms[tpr_to_pmd[i3]]
        a4 = struct.atoms[tpr_to_pmd[i4]]
        if functype == 2 and len(params) >= 2:
            # Harmonic improper (CHARMM-style): psi0 (deg), kpsi (kJ/mol/rad^2)
            psi0, kpsi = float(params[0]), float(params[1])
            d = pmd.Improper(a1, a2, a3, a4)
            d.type = _get_improper_type(psi0, kpsi)
            struct.impropers.append(d)
        elif functype == 4 and len(params) >= 3:
            # Periodic improper: same functional form as a proper dihedral
            phi0, kphi, mult = float(params[0]), float(params[1]), params[2]
            d = pmd.Dihedral(a1, a2, a3, a4, improper=True)
            d.type = _get_dihedral_type(phi0, kphi, mult)
            struct.dihedrals.append(d)

    # --- CMAP (CHARMM backbone grid correction): NOT parsed ---
    # CHARMM-family force fields (used throughout this package's KCX/protein
    # templates) apply a 2D grid correction (CMAP) to the phi/psi backbone
    # dihedral energy. Mapping GROMACS' CMAP grid storage to ParmEd's
    # CmapType/Cmap objects is not implemented here; any CMAP terms present in
    # the TPR are silently absent from the resulting OpenMM System. This is a
    # known, deliberate limitation of Native TPR mode (unlike .prmtop native
    # mode, which is unaffected) and should be treated as an approximation for
    # any energetics that are sensitive to backbone conformational strain.
    try:
        cmaps = reader.get_bonded('cmaps')
        if cmaps is not None and len(cmaps) > 0:
            print(f"[TprLoader] WARNING: {len(cmaps)} CMAP backbone correction term(s) found in TPR "
                  f"but are NOT applied in Native TPR mode (CMAP parsing is unimplemented). "
                  f"Bonded internal energies from this system will omit the CHARMM CMAP correction. "
                  f"For CMAP-accurate results, use .prmtop Native mode instead.")
    except Exception:
        pass

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
    
    # Needs a box so trajectory files don't fail length checks if they have unitcells.
    # GROMACS stores the box as a row-major lower-triangular matrix of vectors
    # v1=(ax,0,0), v2=(bx,by,0), v3=(cx,cy,cz) (nm) -- NOT necessarily
    # orthorhombic. Truncated-octahedron/rhombic-dodecahedron boxes (common
    # choices for protein-ligand production MD to minimize solvent padding)
    # have nonzero off-diagonal components, so the (a,b,c,alpha,beta,gamma)
    # unit cell must be derived from the full vectors, not just the diagonal.
    try:
        box = reader.get_xvf('box')
        if len(box) == 3:
            v1 = np.array(box[0], dtype=float)
            v2 = np.array(box[1], dtype=float)
            v3 = np.array(box[2], dtype=float)
            a = np.linalg.norm(v1)
            b = np.linalg.norm(v2)
            c = np.linalg.norm(v3)
            if a > 0 and b > 0 and c > 0:
                alpha = np.degrees(np.arccos(np.clip(np.dot(v2, v3) / (b * c), -1.0, 1.0)))
                beta = np.degrees(np.arccos(np.clip(np.dot(v1, v3) / (a * c), -1.0, 1.0)))
                gamma = np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (a * b), -1.0, 1.0)))
                struct.box = [a * 10.0, b * 10.0, c * 10.0, alpha, beta, gamma]
    except Exception:
        pass

    return struct
