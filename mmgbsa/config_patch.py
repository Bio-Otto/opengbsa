    def _normalize_config(self):
        """
        Normalize configuration to standard format.
        Handles legacy 'input' section alias to 'input_files'.
        """
        if 'input_files' not in self.config and 'input' in self.config:
            logger.info("Normalizing legacy 'input' section to 'input_files'")
            inp = self.config['input']
            
            # Create input_files section
            self.config['input_files'] = {}
            
            # Map keys
            mapping = {
                'topology': 'complex_pdb',
                'trajectory': 'trajectory',
                'ligand_mol': 'ligand_mol',
                'ligand_pdb': 'ligand_pdb',
                'receptor_topology': 'receptor_topology',
                'ligand_topology': 'ligand_topology',
                'solvated_topology': 'solvated_topology'
            }
            
            for legacy_key, new_key in mapping.items():
                if legacy_key in inp:
                    self.config['input_files'][new_key] = inp[legacy_key]
            
            # Validation requires specific keys, ensure defaults if missing
            # Only if they were optional in legacy but required in new schema
            # But legacy 'input' usually has 'topology' and 'trajectory'.
