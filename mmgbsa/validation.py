"""
Validation module for topology parameters and system consistency checks.

This module provides validation utilities to detect non-standard parameters
and potential issues in topology files, helping users avoid common pitfalls.
"""

import logging
from typing import List, Dict, Optional, Any
import parmed as pmd

logger = logging.getLogger(__name__)


class ValidationWarning:
    """Represents a validation warning or error"""
    
    def __init__(
        self,
        warning_type: str,
        severity: str,
        message: str,
        impact: str,
        suggestion: str = "",
        documentation: str = ""
    ):
        self.type = warning_type
        self.severity = severity  # 'INFO', 'WARNING', 'ERROR'
        self.message = message
        self.impact = impact
        self.suggestion = suggestion
        self.documentation = documentation
    
    def __str__(self):
        lines = [f"{self.severity}: {self.message}"]
        if self.impact:
            lines.append(f"  Impact: {self.impact}")
        if self.suggestion:
            lines.append(f"  Suggestion: {self.suggestion}")
        if self.documentation:
            lines.append(f"  Documentation: {self.documentation}")
        return "\n".join(lines)


class TopologyValidator:
    """Validate topology parameters against AMBER standards"""
    
    # Standard AMBER parameter values
    STANDARD_SCNB = 2.0
    STANDARD_SCEE = 1.2
    
    @staticmethod
    def validate_scnb_scee(structure: pmd.Structure) -> List[ValidationWarning]:
        """
        Check SCNB/SCEE scaling factors against AMBER standards.
        
        Args:
            structure: ParmEd structure to validate
            
        Returns:
            List of validation warnings
        """
        warnings = []
        
        if not hasattr(structure, 'parm_data'):
            return warnings
        
        parm_data = structure.parm_data
        
        # Check SCNB (1-4 VDW scaling)
        if 'SCNB_SCALE_FACTOR' in parm_data:
            scnb_values = parm_data['SCNB_SCALE_FACTOR']
            
            if scnb_values:
                avg_scnb = sum(scnb_values) / len(scnb_values)
                
                # Check if all values are consistent
                if not all(abs(v - scnb_values[0]) < 0.01 for v in scnb_values):
                    warnings.append(ValidationWarning(
                        warning_type='SCNB_INCONSISTENT',
                        severity='WARNING',
                        message='SCNB values are not consistent across dihedrals',
                        impact='May indicate corrupted topology file',
                        suggestion='Regenerate topology file with tleap or check source'
                    ))
                
                # Check if non-standard
                if not all(abs(v - TopologyValidator.STANDARD_SCNB) < 0.01 for v in scnb_values):
                    warnings.append(ValidationWarning(
                        warning_type='SCNB_NON_STANDARD',
                        severity='WARNING',
                        message=f'SCNB={avg_scnb:.2f} (AMBER standard is {TopologyValidator.STANDARD_SCNB})',
                        impact='VDW energies will differ from standard AMBER/MMPBSA.py by ~10-15 kcal/mol',
                        suggestion='Consider regenerating topology with standard AMBER parameters, or accept the difference',
                        documentation='See docs/parameters/force_field_parameters.md for details'
                    ))
        
        # Check SCEE (1-4 Electrostatic scaling)
        if 'SCEE_SCALE_FACTOR' in parm_data:
            scee_values = parm_data['SCEE_SCALE_FACTOR']
            
            if scee_values:
                avg_scee = sum(scee_values) / len(scee_values)
                
                # Check if non-standard
                if not all(abs(v - TopologyValidator.STANDARD_SCEE) < 0.01 for v in scee_values):
                    warnings.append(ValidationWarning(
                        warning_type='SCEE_NON_STANDARD',
                        severity='WARNING',
                        message=f'SCEE={avg_scee:.2f} (AMBER standard is {TopologyValidator.STANDARD_SCEE})',
                        impact='Electrostatic energies may differ from standard AMBER',
                        suggestion='Consider regenerating topology with standard AMBER parameters'
                    ))
        
        return warnings
    
    @staticmethod
    def validate_atom_counts(
        complex_structure: pmd.Structure,
        receptor_structure: pmd.Structure,
        ligand_structure: pmd.Structure
    ) -> List[ValidationWarning]:
        """
        Verify atom count consistency between complex and components.
        
        Args:
            complex_structure: Complex structure
            receptor_structure: Receptor structure
            ligand_structure: Ligand structure
            
        Returns:
            List of validation warnings
        """
        warnings = []
        
        complex_atoms = len(complex_structure.atoms)
        receptor_atoms = len(receptor_structure.atoms)
        ligand_atoms = len(ligand_structure.atoms)
        expected_atoms = receptor_atoms + ligand_atoms
        
        if complex_atoms != expected_atoms:
            warnings.append(ValidationWarning(
                warning_type='ATOM_COUNT_MISMATCH',
                severity='ERROR',
                message=f'Atom count mismatch: complex={complex_atoms}, receptor={receptor_atoms}, ligand={ligand_atoms}, expected={expected_atoms}',
                impact='Calculation will fail or produce incorrect results',
                suggestion='Check structure splitting logic or input files'
            ))
        
        return warnings
    
    @staticmethod
    def validate_force_field(structure: pmd.Structure) -> List[ValidationWarning]:
        """
        Check for common force field issues.
        
        Args:
            structure: ParmEd structure to validate
            
        Returns:
            List of validation warnings
        """
        warnings = []
        
        # Check for missing parameters
        for atom in structure.atoms:
            if atom.epsilon == 0 and atom.rmin == 0:
                warnings.append(ValidationWarning(
                    warning_type='MISSING_VDW_PARAMS',
                    severity='WARNING',
                    message=f'Atom {atom.name} (residue {atom.residue.name}) has zero VDW parameters',
                    impact='May indicate missing force field parameters',
                    suggestion='Check if all residues are properly parameterized'
                ))
                break  # Only report once
        
        return warnings
    
    @staticmethod
    def validate_system_sanity(
        binding_energy: float,
        std_dev: float,
        components: Optional[Dict[str, float]] = None
    ) -> List[ValidationWarning]:
        """
        Check if calculated results are reasonable.
        
        Args:
            binding_energy: Calculated binding energy
            std_dev: Standard deviation
            components: Optional dict of energy components
            
        Returns:
            List of validation warnings
        """
        warnings = []
        
        # Check for unusually large binding energy
        if abs(binding_energy) > 100:
            warnings.append(ValidationWarning(
                warning_type='UNUSUAL_BINDING_ENERGY',
                severity='ERROR',
                message=f'Binding energy ({binding_energy:.2f} kcal/mol) is unusually large',
                impact='Likely indicates calculation error or parameter issue',
                suggestion='Check topology parameters, especially SCNB/SCEE values',
                documentation='See docs/troubleshooting.md#unusual-energies'
            ))
        
        # Check for high variance
        if std_dev > 20:
            warnings.append(ValidationWarning(
                warning_type='HIGH_VARIANCE',
                severity='WARNING',
                message=f'High standard deviation ({std_dev:.2f} kcal/mol) indicates instability',
                impact='Results may not be reliable',
                suggestion='Consider using more frames or checking trajectory stability',
                documentation='See docs/troubleshooting.md#high-variance'
            ))
        
        # Check component sanity if provided
        if components:
            # VDW should be negative
            if components.get('vdw', 0) > 0:
                warnings.append(ValidationWarning(
                    warning_type='POSITIVE_VDW',
                    severity='WARNING',
                    message=f'VDW energy is positive ({components["vdw"]:.2f} kcal/mol)',
                    impact='Unusual, may indicate repulsive interactions or calculation error',
                    suggestion='Check structure for clashes or parameter issues'
                ))
            
            # GB should typically be positive (unfavorable desolvation)
            if components.get('gb', 0) < -50:
                warnings.append(ValidationWarning(
                    warning_type='UNUSUAL_GB',
                    severity='WARNING',
                    message=f'GB energy is very negative ({components["gb"]:.2f} kcal/mol)',
                    impact='May indicate GB model issue',
                    suggestion='Check GB model selection and radii'
                ))
        
        return warnings


class ValidationReport:
    """Collects and reports validation warnings"""
    
    def __init__(self):
        self.warnings: List[ValidationWarning] = []
    
    def add_warnings(self, warnings: List[ValidationWarning]):
        """Add multiple warnings"""
        self.warnings.extend(warnings)
    
    def add_warning(self, warning: ValidationWarning):
        """Add single warning"""
        self.warnings.append(warning)
    
    def has_errors(self) -> bool:
        """Check if any errors exist"""
        return any(w.severity == 'ERROR' for w in self.warnings)
    
    def has_warnings(self) -> bool:
        """Check if any warnings exist"""
        return any(w.severity == 'WARNING' for w in self.warnings)
    
    def print_summary(self):
        """Print formatted validation summary"""
        if not self.warnings:
            logger.info("✅ All validation checks passed")
            return
        
        # Group by severity
        errors = [w for w in self.warnings if w.severity == 'ERROR']
        warnings = [w for w in self.warnings if w.severity == 'WARNING']
        infos = [w for w in self.warnings if w.severity == 'INFO']
        
        # Print errors
        if errors:
            print("\n" + "="*70)
            print("❌ VALIDATION ERRORS")
            print("="*70)
            for err in errors:
                print(f"\n{err}")
        
        # Print warnings
        if warnings:
            print("\n" + "="*70)
            print("⚠️  VALIDATION WARNINGS")
            print("="*70)
            for warn in warnings:
                print(f"\n{warn}")
        
        # Print infos
        if infos:
            print("\n" + "="*70)
            print("ℹ️  VALIDATION INFO")
            print("="*70)
            for info in infos:
                print(f"\n{info}")
        
        print("\n" + "="*70)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'total_warnings': len(self.warnings),
            'errors': len([w for w in self.warnings if w.severity == 'ERROR']),
            'warnings': len([w for w in self.warnings if w.severity == 'WARNING']),
            'infos': len([w for w in self.warnings if w.severity == 'INFO']),
            'details': [
                {
                    'type': w.type,
                    'severity': w.severity,
                    'message': w.message,
                    'impact': w.impact,
                    'suggestion': w.suggestion
                }
                for w in self.warnings
            ]
        }
