    def create_combined_system(self, protein_system, ligand_system):
        """
        Create a combined OpenMM system containing both protein and ligand
        for accurate MM/GBSA interaction energy calculation.
        
        This method merges two separately parametrized systems (protein and ligand)
        into a single system that can calculate the total energy including
        protein-ligand interactions.
        
        Parameters:
        -----------
        protein_system : openmm.System
            Parametrized protein system with GB forces
        ligand_system : openmm.System
            Parametrized ligand system with GB forces
            
        Returns:
        --------
        openmm.System : Combined system with protein + ligand
        """
        import copy
        
        print("Creating combined protein+ligand system for MM/GBSA...")
        
        # Create new empty system
        combined_system = openmm.System()
        
        # Add all particles from protein
        n_protein = protein_system.getNumParticles()
        for i in range(n_protein):
            mass = protein_system.getParticleMass(i)
            combined_system.addParticle(mass)
        
        # Add all particles from ligand
        n_ligand = ligand_system.getNumParticles()
        for i in range(n_ligand):
            mass = ligand_system.getParticleMass(i)
            combined_system.addParticle(mass)
        
        print(f"  Combined system: {n_protein} protein + {n_ligand} ligand = {combined_system.getNumParticles()} total particles")
        
        # Offsets for indexing
        protein_offset = 0
        ligand_offset = n_protein
        
        # Copy forces from both systems
        # We need to handle each force type separately because they have different APIs
        
        protein_forces = list(protein_system.getForces())
        ligand_forces = list(ligand_system.getForces())
        
        print(f"  Merging {len(protein_forces)} protein forces and {len(ligand_forces)} ligand forces...")
        
        # Strategy: For GB-based MM/GBSA, we need to merge the GB forces
        # The key forces are:
        # 1. CustomGBForce (GB solvation)
        # 2. CustomNonbondedForce (screening/interactions)
        # 3. CustomBondForce (surface area)
        
        # Find and merge GB forces
        protein_gb_force = None
        ligand_gb_force = None
        protein_screening_force = None
        ligand_screening_force = None
        protein_sa_force = None
        ligand_sa_force = None
        
        for force in protein_forces:
            if isinstance(force, openmm.CustomGBForce):
                protein_gb_force = force
            elif isinstance(force, openmm.CustomNonbondedForce):
                protein_screening_force = force
            elif isinstance(force, openmm.CustomBondForce):
                protein_sa_force = force
        
        for force in ligand_forces:
            if isinstance(force, openmm.CustomGBForce):
                ligand_gb_force = force
            elif isinstance(force, openmm.CustomNonbondedForce):
                ligand_screening_force = force
            elif isinstance(force, openmm.CustomBondForce):
                ligand_sa_force = force
        
        # Merge CustomGBForce (GB solvation)
        if protein_gb_force and ligand_gb_force:
            combined_gb = self._merge_custom_gb_forces(protein_gb_force, ligand_gb_force, 
                                                       protein_offset, ligand_offset)
            combined_system.addForce(combined_gb)
            print("  ✓ Merged CustomGBForce")
        
        # Merge CustomNonbondedForce (screening/interactions)
        if protein_screening_force and ligand_screening_force:
            combined_screening = self._merge_custom_nonbonded_forces(
                protein_screening_force, ligand_screening_force,
                protein_offset, ligand_offset)
            combined_system.addForce(combined_screening)
            print("  ✓ Merged CustomNonbondedForce")
        
        # Merge CustomBondForce (surface area)
        if protein_sa_force and ligand_sa_force:
            combined_sa = self._merge_custom_bond_forces(protein_sa_force, ligand_sa_force,
                                                         protein_offset, ligand_offset)
            combined_system.addForce(combined_sa)
            print("  ✓ Merged CustomBondForce")
        
        print(f"✓ Combined system created with {combined_system.getNumForces()} forces")
        
        return combined_system
    
    def _merge_custom_gb_forces(self, protein_force, ligand_force, protein_offset, ligand_offset):
        """Merge two CustomGBForce objects"""
        # Create new GB force with same parameters as protein force
        combined = openmm.CustomGBForce()
        
        # Copy energy function and parameters from protein force
        # (assuming both have same GB model - OBC2)
        combined.addPerParticleParameter("q")
        combined.addPerParticleParameter("radius")
        combined.addPerParticleParameter("scale")
        combined.addGlobalParameter("solventDielectric", 78.5)
        combined.addGlobalParameter("soluteDielectric", 1.0)
        
        # Copy computed values
        for i in range(protein_force.getNumComputedValues()):
            name, expression, cv_type = protein_force.getComputedValueParameters(i)
            combined.addComputedValue(name, expression, cv_type)
        
        # Copy energy terms
        for i in range(protein_force.getNumEnergyTerms()):
            name, expression, et_type = protein_force.getEnergyTermParameters(i)
            combined.addEnergyTerm(name, expression, et_type)
        
        # Add protein particles
        for i in range(protein_force.getNumParticles()):
            params = protein_force.getParticleParameters(i)
            combined.addParticle(params)
        
        # Add ligand particles
        for i in range(ligand_force.getNumParticles()):
            params = ligand_force.getParticleParameters(i)
            combined.addParticle(params)
        
        combined.setNonbondedMethod(protein_force.getNonbondedMethod())
        combined.setCutoffDistance(protein_force.getCutoffDistance())
        
        return combined
    
    def _merge_custom_nonbonded_forces(self, protein_force, ligand_force, protein_offset, ligand_offset):
        """Merge two CustomNonbondedForce objects"""
        # Create new force with same energy function
        energy_function = protein_force.getEnergyFunction()
        combined = openmm.CustomNonbondedForce(energy_function)
        
        # Copy parameters
        for i in range(protein_force.getNumPerParticleParameters()):
            name = protein_force.getPerParticleParameterName(i)
            combined.addPerParticleParameter(name)
        
        for i in range(protein_force.getNumGlobalParameters()):
            name = protein_force.getGlobalParameterName(i)
            value = protein_force.getGlobalParameterDefaultValue(i)
            combined.addGlobalParameter(name, value)
        
        # Add protein particles
        for i in range(protein_force.getNumParticles()):
            params = protein_force.getParticleParameters(i)
            combined.addParticle(params)
        
        # Add ligand particles
        for i in range(ligand_force.getNumParticles()):
            params = ligand_force.getParticleParameters(i)
            combined.addParticle(params)
        
        # Copy protein exclusions
        for i in range(protein_force.getNumExclusions()):
            p1, p2 = protein_force.getExclusionParticles(i)
            combined.addExclusion(p1 + protein_offset, p2 + protein_offset)
        
        # Copy ligand exclusions with offset
        for i in range(ligand_force.getNumExclusions()):
            p1, p2 = ligand_force.getExclusionParticles(i)
            combined.addExclusion(p1 + ligand_offset, p2 + ligand_offset)
        
        combined.setNonbondedMethod(protein_force.getNonbondedMethod())
        combined.setCutoffDistance(protein_force.getCutoffDistance())
        
        return combined
    
    def _merge_custom_bond_forces(self, protein_force, ligand_force, protein_offset, ligand_offset):
        """Merge two CustomBondForce objects"""
        # Create new force with same energy function
        energy_function = protein_force.getEnergyFunction()
        combined = openmm.CustomBondForce(energy_function)
        
        # Copy parameters
        for i in range(protein_force.getNumPerBondParameters()):
            name = protein_force.getPerBondParameterName(i)
            combined.addPerBondParameter(name)
        
        for i in range(protein_force.getNumGlobalParameters()):
            name = protein_force.getGlobalParameterName(i)
            value = protein_force.getGlobalParameterDefaultValue(i)
            combined.addGlobalParameter(name, value)
        
        # Add protein bonds
        for i in range(protein_force.getNumBonds()):
            p1, p2, params = protein_force.getBondParameters(i)
            combined.addBond(p1 + protein_offset, p2 + protein_offset, params)
        
        # Add ligand bonds with offset
        for i in range(ligand_force.getNumBonds()):
            p1, p2, params = ligand_force.getBondParameters(i)
            combined.addBond(p1 + ligand_offset, p2 + ligand_offset, params)
        
        return combined
