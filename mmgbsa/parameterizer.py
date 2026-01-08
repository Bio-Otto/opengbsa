import logging
import os
import shutil
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
try:
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField as OpenFFForceField
    HAS_OPENFF = True
except ImportError:
    HAS_OPENFF = False

from openmm import app, unit, System
import xml.etree.ElementTree as ET
from xml.dom import minidom

log = logging.getLogger(__name__)

class ResidueParameterizer:
    def __init__(self, work_dir="parameterization_work"):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def extract_and_cap_residue(self, pdb_path, res_name, output_sdf):
        print(f"DEBUG: Extracting {res_name} from {pdb_path}")
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
        if not mol:
            raise ValueError(f"Could not load PDB: {pdb_path}")

        # 1. Identify Residue Atoms
        res_atoms = []
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == res_name:
                res_atoms.append(atom.GetIdx())
        
        if not res_atoms:
            raise ValueError(f"Residue {res_name} not found in {pdb_path}")

        # 1b. Ensure Internal Connectivity (Distance-based)
        conf = mol.GetConformer()
        added_bonds = 0
        em_bond = Chem.EditableMol(mol)
        
        for i in range(len(res_atoms)):
            idx1 = res_atoms[i]
            atom1 = mol.GetAtomWithIdx(idx1)
            if atom1.GetAtomicNum() == 1: continue 
            
            pos1 = conf.GetAtomPosition(idx1)
            for j in range(i + 1, len(res_atoms)):
                idx2 = res_atoms[j]
                atom2 = mol.GetAtomWithIdx(idx2)
                if atom2.GetAtomicNum() == 1: continue 
                
                pos2 = conf.GetAtomPosition(idx2)
                dist = pos1.Distance(pos2)
                
                if dist < 1.9:
                    bond = mol.GetBondBetweenAtoms(idx1, idx2)
                    if not bond:
                        print(f"DEBUG: Adding bond {idx1}-{idx2} dist={dist}")
                        em_bond.AddBond(idx1, idx2, Chem.BondType.SINGLE)
                        added_bonds += 1
                        
        if added_bonds > 0:
            mol = em_bond.GetMol()
            print(f"DEBUG: Inferred {added_bonds} internal bonds.")
        
        # 2. Extract Substructure
        em = Chem.EditableMol(mol)
        all_indices = set(range(mol.GetNumAtoms()))
        to_delete = sorted(list(all_indices - set(res_atoms)), reverse=True)
        for idx in to_delete:
            em.RemoveAtom(idx)
        res_mol = em.GetMol()
        
        # 3. Identify Backbone N and C and Capture Names
        n_idx = -1
        c_idx = -1
        mapping_dict = {} 
        
        for atom in res_mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            name = info.GetName().strip()
            mapping_dict[atom.GetIdx()] = name
            
            if name == "N":
                n_idx = atom.GetIdx()
            elif name == "C":
                c_idx = atom.GetIdx()
        
        # 4. Cap with ACE and NME
        em_cap = Chem.EditableMol(res_mol)
        def add_cap_atom(atomic_num):
            return em_cap.AddAtom(Chem.Atom(atomic_num))

        # ACE
        ace_c_idx = add_cap_atom(6) 
        ace_o_idx = add_cap_atom(8)
        ace_ch3_idx = add_cap_atom(6)
        if n_idx != -1:
            em_cap.AddBond(n_idx, ace_c_idx, Chem.BondType.SINGLE)
        em_cap.AddBond(ace_c_idx, ace_o_idx, Chem.BondType.DOUBLE)
        em_cap.AddBond(ace_c_idx, ace_ch3_idx, Chem.BondType.SINGLE)
        
        # NME
        nme_n_idx = add_cap_atom(7)
        nme_ch3_idx = add_cap_atom(6)
        if c_idx != -1:
            em_cap.AddBond(c_idx, nme_n_idx, Chem.BondType.SINGLE)
        em_cap.AddBond(nme_n_idx, nme_ch3_idx, Chem.BondType.SINGLE)
        
        capped_mol = em_cap.GetMol()
        
        num_res_atoms = res_mol.GetNumAtoms() # Original (Heavy+H from PDB)
        num_capped_atoms = capped_mol.GetNumAtoms() # Original + Caps (Heavy)
        
        
        # CRITICAL FIX: Only add hydrogens to CAP atoms, preserve residue atoms from PDB
        # Strategy: Manually saturate cap atoms with correct number of hydrogens
        em_final = Chem.EditableMol(capped_mol)
        
        # ACE cap needs 3 H on CH3 and none on C=O
        for _ in range(3):
            h_idx = em_final.AddAtom(Chem.Atom(1))  # Add H
            em_final.AddBond(ace_ch3_idx, h_idx, Chem.BondType.SINGLE)
        
        # NME cap needs 1 H on N and 3 H on CH3
        h_idx = em_final.AddAtom(Chem.Atom(1))  # H on N
        em_final.AddBond(nme_n_idx, h_idx, Chem.BondType.SINGLE)
        for _ in range(3):
            h_idx = em_final.AddAtom(Chem.Atom(1))  # Add H on CH3
            em_final.AddBond(nme_ch3_idx, h_idx, Chem.BondType.SINGLE)
        
        mapped_mol = em_final.GetMol()
        
        # Atom name mapping
        final_mapping = {}
        
        # Indices:
        # 0 .. num_res_atoms-1 : Original Core Atoms (Heavy + PDB Hs)
        # num_res_atoms .. num_capped_atoms-1 : Cap Atoms (Heavy)
        # num_capped_atoms .. end : NEW Cap Hydrogens (manually added)
        
        existing_names = set()
        for i in range(num_res_atoms):
            if i in mapping_dict:
                existing_names.add(mapping_dict[i])
        
        h_counter = 1
        
        for atom in mapped_mol.GetAtoms():
            idx = atom.GetIdx()
            
            # 1. Original Core Atoms
            if idx < num_res_atoms:
                atom.SetBoolProp("is_cap", False)
                if idx in mapping_dict:
                     name = mapping_dict[idx]
                     final_mapping[idx] = name
                     atom.SetProp("atom_name", name)
                else:
                     # Should not happen
                     name = f"X{idx}"
                     final_mapping[idx] = name
                     atom.SetProp("atom_name", name)
            
            # 2. Cap Atoms
            elif idx < num_capped_atoms:
                atom.SetBoolProp("is_cap", True)
                atom.SetProp("name", "CAP")
            
            # 3. New Hydrogens
            else:
                # Check neighbor to determine if Core or Cap
                is_core_h = False
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx < num_res_atoms:
                        is_core_h = True
                        break
                
                if is_core_h:
                    atom.SetBoolProp("is_cap", False)
                    # Generate Short Unique Name
                    while True:
                        new_name = f"H{h_counter:02d}" # H01, H02...
                        h_counter += 1
                        if new_name not in existing_names:
                            existing_names.add(new_name)
                            break
                    
                    final_mapping[idx] = new_name
                    atom.SetProp("atom_name", new_name)
                else:
                    atom.SetBoolProp("is_cap", True)
                    atom.SetProp("name", "CAP")

        # 5. Sanitize
        try:
            Chem.SanitizeMol(capped_mol)
        except Exception as e:
            print(f"Sanitization warning: {e}")
        
        # 6. Embed
        AllChem.EmbedMolecule(capped_mol, AllChem.ETKDG())
        
        # Save SDF
        w = Chem.SDWriter(str(output_sdf))
        w.SetForceV3000(True)
        w.write(capped_mol)
        w.close()
        
        # SAVE MAPPING JSON
        json_path = str(output_sdf).replace(".sdf", "_map.json")
        with open(json_path, "w") as f:
            json.dump(final_mapping, f, indent=2)
        
        # 7. Create Protonated Residue PDB (for PDB Surgery)
        protonated_mol = Chem.RWMol(mapped_mol)
        
        atoms_to_remove = []
        for atom in protonated_mol.GetAtoms():
            if atom.HasProp("is_cap") and atom.GetBoolProp("is_cap"):
                atoms_to_remove.append(atom.GetIdx())
            else:
                # Ensure PDB Info has correct Name (for XML matching)
                name = atom.GetProp("atom_name") if atom.HasProp("atom_name") else atom.GetSymbol()
                info = atom.GetPDBResidueInfo()
                if not info:
                    info = Chem.AtomPDBResidueInfo()
                    info.SetResidueName(res_name)
                    info.SetIsHeteroAtom(False)
                    atom.SetPDBResidueInfo(info)
                
                # Align Name for PDB (4 chars)
                # If name is 1-3 chars (e.g. "N"), PDB usually writes " N  " (padded 1 left, left aligned)
                # If 4 chars, writes "XXXX"
                # RDKit SetName writes logic?
                # RDKit PDBWriter tries to be smart.
                # Just set the name.
                info.SetName(name)
                atom.SetPDBResidueInfo(info)

        # Remove Caps in reverse order
        atoms_to_remove.sort(reverse=True)
        for idx in atoms_to_remove:
            protonated_mol.RemoveAtom(idx)

        protonated_pdb_path = str(output_sdf).replace("_capped.sdf", "_protonated.pdb")
        Chem.MolToPDBFile(protonated_mol, protonated_pdb_path)
        print(f"DEBUG: Saved protonated residue PDB to {protonated_pdb_path}")
        
        return output_sdf, protonated_pdb_path

    def parameterize_with_openff(self, sdf_path, output_xml_path, res_name, charge_method=None):
        if not HAS_OPENFF:
            raise ImportError("OpenFF Toolkit not installed.")
            
        print(f"DEBUG: Parameterizing {res_name}...")
        molecule = Molecule.from_file(str(sdf_path))
        print(f"DEBUG: OpenFF Mol has {molecule.n_atoms} atoms and {len(molecule.bonds)} bonds.")
        
        charge_kwargs = {}
        if charge_method:
             print(f"DEBUG: Assigning charges using {charge_method}")
             molecule.assign_partial_charges(partial_charge_method=charge_method)
             charge_kwargs = {'charge_from_molecules': [molecule]}
        
        ff = OpenFFForceField("openff-2.0.0.offxml")
        openmm_system = ff.create_openmm_system(molecule.to_topology(), **charge_kwargs)
        
        # Load JSON map
        json_path = str(sdf_path).replace(".sdf", "_map.json")
        atom_map = {}
        if os.path.exists(json_path):
            with open(json_path, "r") as f:
                atom_map = json.load(f)
                atom_map = {int(k): v for k, v in atom_map.items()}
        
        self._convert_system_to_residue_xml(openmm_system, molecule, res_name, output_xml_path, atom_map, str(sdf_path))
        return output_xml_path

    def _convert_system_to_residue_xml(self, system, openff_mol, res_name, output_path, atom_map=None, sdf_path=None):
        # Reload SDF to ensure index matching if needed
        rdmol = openff_mol.to_rdkit()
        if sdf_path:
             ref = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
             if ref and ref.GetNumAtoms() == openff_mol.n_atoms:
                 rdmol = ref

        core_indices = []
        for i, atom in enumerate(rdmol.GetAtoms()):
            if atom_map and i in atom_map:
                core_indices.append(i)
            elif atom.HasProp("is_cap") and atom.GetBoolProp("is_cap"):
                 continue
            if not atom_map and not (atom.HasProp("is_cap") and atom.GetBoolProp("is_cap")):
                 core_indices.append(i)
            
        core_set = set(core_indices)
        print(f"DEBUG: Core Set size: {len(core_set)}")

        root = ET.Element("ForceField")
        atom_types = ET.SubElement(root, "AtomTypes")
        residues = ET.SubElement(root, "Residues")
        res_elem = ET.SubElement(residues, "Residue", name=res_name)
        
        hb_force_elem = ET.SubElement(root, "HarmonicBondForce")
        ha_force_elem = ET.SubElement(root, "HarmonicAngleForce")
        pt_force_elem = ET.SubElement(root, "PeriodicTorsionForce")
        nb_force_elem = ET.SubElement(root, "NonbondedForce", coulomb14scale="0.833333333333", lj14scale="0.5") 
        
        sys_forces = {f.__class__.__name__: f for f in system.getForces()}
        nb_force = sys_forces['NonbondedForce']
        
        atom_index_to_type = {}
        atom_index_to_name = {}
        sys_idx_to_xml_idx = {} 
        
        xml_counter = 0
        for i in core_indices:
            atom = rdmol.GetAtomWithIdx(i)
            
            if atom_map and i in atom_map:
                atom_name = atom_map[i]
            else:
                info = atom.GetPDBResidueInfo()
                atom_name = info.GetName().strip() if info else f"A{i}"
            
            atom_index_to_name[i] = atom_name
            sys_idx_to_xml_idx[i] = xml_counter
            xml_counter += 1
            
            type_name = f"{res_name}-{atom_name}"
            atom_index_to_type[i] = type_name
            
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            
            at_elem = ET.SubElement(atom_types, "Type", 
                          name=type_name, 
                          element=atom.GetSymbol(),
                          mass=str(atom.GetMass()))
            at_elem.set("class", type_name) 
            
            ET.SubElement(nb_force_elem, "Atom",
                          type=type_name,
                          charge=str(charge.value_in_unit(unit.elementary_charge)),
                          sigma=str(sigma.value_in_unit(unit.nanometers)),
                          epsilon=str(epsilon.value_in_unit(unit.kilojoules_per_mole)))
            
            ET.SubElement(res_elem, "Atom", name=atom_name, type=type_name)

        # Bonds
        count_bonds = 0
        print(f"DEBUG: Checking {len(openff_mol.bonds)} OpenFF bonds against Core Set")
        for bond in openff_mol.bonds:
            i = bond.atom1_index
            j = bond.atom2_index
            
            if i in core_set and j in core_set:
                b_elem = ET.SubElement(res_elem, "Bond")
                b_elem.set("from", str(sys_idx_to_xml_idx[i])) 
                b_elem.set("to", str(sys_idx_to_xml_idx[j]))  
                count_bonds += 1
            
            elif (i in core_set and j not in core_set) or (j in core_set and i not in core_set):
                core_idx = i if i in core_set else j
                name = atom_index_to_name[core_idx]
                ET.SubElement(res_elem, "ExternalBond", atomName=name)
        
        print(f"DEBUG: Added {count_bonds} internal bonds to XML")

        # Forces
        hb_force = sys_forces.get('HarmonicBondForce')
        if hb_force:
            for k in range(hb_force.getNumBonds()):
                idx1, idx2, length, k_const = hb_force.getBondParameters(k)
                if idx1 in core_set and idx2 in core_set:
                    t1 = atom_index_to_type[idx1]
                    t2 = atom_index_to_type[idx2]
                    ET.SubElement(hb_force_elem, "Bond", 
                                  class1=t1, class2=t2, 
                                  length=str(length.value_in_unit(unit.nanometers)),
                                  k=str(k_const.value_in_unit(unit.kilojoules_per_mole/unit.nanometer**2)))

        ha_force = sys_forces.get('HarmonicAngleForce')
        if ha_force:
            for k in range(ha_force.getNumAngles()):
                idx1, idx2, idx3, angle, k_const = ha_force.getAngleParameters(k)
                if idx1 in core_set and idx2 in core_set and idx3 in core_set:
                    ET.SubElement(ha_force_elem, "Angle", 
                                  class1=atom_index_to_type[idx1], 
                                  class2=atom_index_to_type[idx2], 
                                  class3=atom_index_to_type[idx3],
                                  angle=str(angle.value_in_unit(unit.radians)),
                                  k=str(k_const.value_in_unit(unit.kilojoules_per_mole/unit.radian**2)))

        pt_force = sys_forces.get('PeriodicTorsionForce')
        if pt_force:
            for k in range(pt_force.getNumTorsions()):
                idx1, idx2, idx3, idx4, period, phase, k_const = pt_force.getTorsionParameters(k)
                if {idx1, idx2, idx3, idx4}.issubset(core_set):
                     ET.SubElement(pt_force_elem, "Proper", 
                                  class1=atom_index_to_type[idx1], 
                                  class2=atom_index_to_type[idx2], 
                                  class3=atom_index_to_type[idx3],
                                  class4=atom_index_to_type[idx4],
                                  periodicity1=str(period),
                                  phase1=str(phase.value_in_unit(unit.radians)),
                                  k1=str(k_const.value_in_unit(unit.kilojoules_per_mole)))

        xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        with open(output_path, "w") as f:
            f.write(xml_str)
        log.info(f"Generated XML at {output_path}")

    def simple_run(self, pdb_path, res_name, output_xml, charge_method=None):
        sdf = self.work_dir / f"{res_name}_capped.sdf"
        _, protonated_pdb = self.extract_and_cap_residue(pdb_path, res_name, sdf)
        self.parameterize_with_openff(sdf, output_xml, res_name, charge_method=charge_method)
        return output_xml, protonated_pdb
