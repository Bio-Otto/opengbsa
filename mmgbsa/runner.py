#!/usr/bin/env python3
"""
MM/GBSA Analysis Runner with YAML Configuration
Single entry point for running complete MM/GBSA analysis with all features
"""

import yaml
import argparse
import sys
from pathlib import Path
import time
from datetime import datetime
import warnings

# Suppress specific warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="Unable to load toolkit 'OpenEye Toolkit'")
warnings.filterwarnings("ignore", message="importing 'simtk.openmm' is deprecated")

# Import your existing modules
from .core import GBSACalculator
from .entropy import run_ultra_robust_nma
from .decomposition import PerResidueDecomposition
from .logger import ToolLogger

# Initialize global logger
log = ToolLogger()

class MMGBSARunner:
    """
    Main runner class for MM/GBSA analysis with YAML configuration
    """
    
    def __init__(self, config_file, output_dir=None):
        """
        Initialize runner with YAML configuration
        
        Parameters:
        -----------
        config_file : str or dict
            Path to YAML configuration file or config dictionary
        output_dir : str, optional
            Output directory path
        """
        if isinstance(config_file, str):
            self.config_file = config_file
            self.config = self._load_config()
        else:
            self.config_file = None
            self.config = config_file
        
        self.results = {}
        self.output_dir = output_dir
        
    def _load_config(self):
        """Load and validate YAML configuration"""
        try:
            with open(self.config_file, 'r') as f:
                config = yaml.safe_load(f)
            
            # Validate required sections
            required_sections = ['input_files', 'analysis_settings']
            for section in required_sections:
                if section not in config:
                    raise ValueError(f"Missing required section: {section}")
            
            # Set defaults for optional sections
            if 'output_settings' not in config:
                config['output_settings'] = {}
            if 'advanced_settings' not in config:
                config['advanced_settings'] = {}
            
            return config
            
        except Exception as e:
            log.error(f"Error loading configuration: {e}")
            sys.exit(1)
    
    def _validate_input_files(self):
        """Validate all input files exist"""
        input_files = self.config['input_files']
        missing_files = []
        
        for file_type, file_path in input_files.items():
            # Skip None/null values
            if file_path is None:
                continue
            
            if not Path(file_path).exists():
                missing_files.append(f"{file_type}: {file_path}")
        
        if missing_files:
            log.error("Missing input files:")
            for missing in missing_files:
                log.error(f"  • {missing}")
            return False
        
        log.success("All required input files found")
        return True
    
    def _select_frames(self, trajectory_length):
        """
        Select frames based on configuration parameters
        """
        analysis_settings = self.config['analysis_settings']
        
        # Get frame selection parameters
        max_frames = analysis_settings.get('max_frames', None)
        frame_start = analysis_settings.get('frame_start', None)
        frame_end = analysis_settings.get('frame_end', None)
        frame_stride = analysis_settings.get('frame_stride', None)
        frame_selection = analysis_settings.get('frame_selection', 'sequential')
        random_seed = analysis_settings.get('random_seed', 42)
        
        # Determine frame range
        if frame_start is None:
            frame_start = 0
        if frame_end is None:
            frame_end = trajectory_length
        
        # Validate frame range
        frame_start = max(0, min(frame_start, trajectory_length - 1))
        frame_end = max(frame_start + 1, min(frame_end, trajectory_length))
        
        log.info("Frame selection parameters:")
        log.result("Trajectory length", trajectory_length)
        log.result("Frame range", f"{frame_start} to {frame_end}")
        log.result("Selection method", frame_selection)
        
        # ... (Selection logic remains same) ...
        # Generate frame indices based on selection method
        if frame_selection == "sequential":
            if frame_stride is None or frame_stride <= 0:
                frame_indices = list(range(frame_start, frame_end))
            else:
                frame_indices = list(range(frame_start, frame_end, frame_stride))
                
        elif frame_selection == "equidistant":
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            step = (frame_end - frame_start) / max_frames
            frame_indices = [frame_start + int(i * step) for i in range(max_frames)]
            frame_indices = [min(idx, frame_end - 1) for idx in frame_indices]
            
        elif frame_selection == "random":
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            import random
            random.seed(random_seed)
            available_frames = list(range(frame_start, frame_end))
            if frame_stride is not None and frame_stride > 1:
                available_frames = available_frames[::frame_stride]
            
            frame_indices = random.sample(available_frames, min(max_frames, len(available_frames)))
            frame_indices.sort()  # Keep chronological order
            
        else:
            raise ValueError(f"Unknown frame selection method: {frame_selection}")
        
        # Apply max_frames limit
        if max_frames is not None and max_frames > 0:
            frame_indices = frame_indices[:max_frames]
        
        # Remove duplicates and sort
        frame_indices = sorted(list(set(frame_indices)))
        
        log.success(f"Selected {len(frame_indices)} frames for analysis")
        
        return frame_indices
    
    def _create_output_directory(self):
        """Create output directory with timestamp"""
        if self.output_dir is None:
            output_dir = self.config['output_settings'].get('output_directory', 'mmgbsa_results')
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = Path(output_dir) / f"analysis_{timestamp}"
        else:
            output_path = Path(self.output_dir)
        
        output_path.mkdir(parents=True, exist_ok=True)
        
        log.info(f"Output directory established: {output_path}")
        return output_path
    
    def _setup_mmgbsa_calculator(self):
        """Setup MM/GBSA calculator with configuration"""
        analysis_settings = self.config['analysis_settings']
        forcefield_settings = self.config.get('forcefield_settings', {})
        
        calculator = GBSACalculator(
            temperature=analysis_settings.get('temperature', 300),
            verbose=analysis_settings.get('verbose', 1),
            gb_model=analysis_settings.get('gb_model', 'OBC2'),
            salt_concentration=analysis_settings.get('salt_concentration', 0.15),
            use_cache=analysis_settings.get('use_cache', True),
            parallel_processing=analysis_settings.get('parallel_processing', False),
            max_workers=analysis_settings.get('max_workers', None),
            protein_forcefield=forcefield_settings.get('protein_forcefield', 'amber'),
            charge_method=analysis_settings.get('charge_method', 'am1bcc'),
            solute_dielectric=analysis_settings.get('solute_dielectric', 1.0),
            solvent_dielectric=analysis_settings.get('solvent_dielectric', 78.5),
            entropy_method=analysis_settings.get('entropy_method', 'none'),
            decomposition_method=analysis_settings.get('decomposition_method', 'full')
        )
        
        return calculator
    
    def run_analysis(self):
        """Run complete MM/GBSA analysis pipeline"""
        
        log.header("MM/GBSA ANALYSIS RUNNER")
        
        start_time = time.time()
        
        # Step 1: Validate configuration
        log.section("STEP 1: Configuration Validation")
        
        if not self._validate_input_files():
            return None
        
        # Step 2: Create output directory
        log.section("STEP 2: Output Setup")
        
        output_dir = self._create_output_directory()
        
        # Step 3: Setup calculator
        log.section("STEP 3: Calculator Setup")
        
        calculator = self._setup_mmgbsa_calculator()
        
        # Step 4: Run MM/GBSA analysis
        log.section("STEP 4: MM/GBSA Analysis")
        
        input_files = self.config['input_files']
        analysis_settings = self.config['analysis_settings']
        
        # Get frame selection parameters for reporting
        frame_params = {
            'max_frames': analysis_settings.get('max_frames', 50),
            'frame_start': analysis_settings.get('frame_start', None),
            'frame_end': analysis_settings.get('frame_end', None),
            'frame_stride': analysis_settings.get('frame_stride', None),
            'frame_selection': analysis_settings.get('frame_selection', 'sequential')
        }
        
        mmgbsa_results = calculator.run_enhanced(
            ligand_mol=input_files.get('ligand_mol'),
            complex_pdb=input_files['complex_pdb'],
            xtc_file=input_files['trajectory'],
            ligand_pdb=input_files.get('ligand_pdb'),
            max_frames=analysis_settings.get('max_frames', 50),
            energy_decomposition=analysis_settings.get('energy_decomposition', False),
            qha_analyze_complex=analysis_settings.get('qha_analyze_complex', False),
            output_dir=output_dir,
            ligand_selection=analysis_settings.get('ligand_selection'),
            receptor_selection=analysis_settings.get('receptor_selection')
        )
        
        if not mmgbsa_results:
            log.error("MM/GBSA analysis failed!")
            return None
        
        self.results['mmgbsa'] = mmgbsa_results
        log.success(f"MM/GBSA (Standard) completed: {mmgbsa_results['mean_binding_energy']:.2f} ± {mmgbsa_results['std_error']:.2f} kcal/mol")

        if mmgbsa_results.get('parameterized_residues'):
            log.result("Parameterized Residues", ", ".join(mmgbsa_results['parameterized_residues']))

        # DUAL MODE: Dimer Interface Analysis
        if analysis_settings.get('dimer_mode', False):
            log.section("STEP 4.5: Dimer Interface Analysis (Dual Mode)")
            log.info("Calculating interface stability between defined subunits...")
            
            receptor_subunits = analysis_settings.get('receptor_subunits', {})
            unit_a = receptor_subunits.get('unit_a')
            unit_b = receptor_subunits.get('unit_b')
            
            if unit_a and unit_b:
                try:
                    interface_dir = output_dir / 'interface_stability'
                    interface_dir.mkdir(exist_ok=True)
                    
                    log.info(f"  • Receptor Mask: {unit_a}")
                    log.info(f"  • Ligand Mask:   {unit_b} (Treated as ligand for calculation)")
                    
                    interface_results = calculator.run_enhanced(
                        ligand_mol=None, # Force Protein-Protein mode
                        complex_pdb=input_files['complex_pdb'],
                        xtc_file=input_files['trajectory'],
                        ligand_pdb=None,
                        max_frames=analysis_settings.get('max_frames', 50),
                        energy_decomposition=False, # Keep it simple for interface
                        qha_analyze_complex=False,
                        output_dir=interface_dir,
                        ligand_selection=unit_b,
                        receptor_selection=unit_a
                    )
                    
                    if interface_results:
                        self.results['interface_stability'] = interface_results
                        log.success(f"Interface Stability: {interface_results['mean_binding_energy']:.2f} ± {interface_results['std_error']:.2f} kcal/mol")
                        
                        # Compare energies
                        binding_E = mmgbsa_results['mean_binding_energy']
                        interface_E = interface_results['mean_binding_energy']
                        log.result("Ligand Binding Energy", f"{binding_E:.2f}", "kcal/mol")
                        log.result("Dimer Interface Energy", f"{interface_E:.2f}", "kcal/mol")
                except Exception as e:
                    log.error(f"Dimer Interface Analysis Failed: {e}")
                    log.warning("Continuing with rest of analysis...")
            else:
                log.warning("dimer_mode is enabled but 'receptor_subunits' (unit_a, unit_b) are missing in config!")
            
        # Step 5: Run entropy analysis (if enabled)
        if analysis_settings.get('run_entropy_analysis', False):
            log.section("STEP 5: Entropy Analysis")
            
            entropy_results = self._run_entropy_analysis(calculator, input_files)
            if entropy_results:
                self.results['entropy'] = entropy_results
                log.success("Entropy analysis completed")
        
        # Step 6: Run per-residue decomposition (if enabled)
        if analysis_settings.get('run_per_residue_decomposition', False) or analysis_settings.get('per_residue_decomposition', False):
            log.section("STEP 6: Per-Residue Decomposition")
            
            # Get frame parameters for per-residue decomposition
            frame_params = {
                'frame_start': analysis_settings.get('frame_start', None),
                'frame_end': analysis_settings.get('frame_end', None),
                'frame_stride': analysis_settings.get('frame_stride', None),
                'frame_selection': analysis_settings.get('frame_selection', 'sequential'),
                'random_seed': analysis_settings.get('random_seed', 42)
            }
            
            decomp_results = self._run_per_residue_decomposition(
                calculator, input_files, analysis_settings, frame_params, output_dir
            )
            if decomp_results:
                self.results['decomposition'] = decomp_results
                log.success("Per-residue decomposition completed")
        
        # Step 7: Generate final report
        log.section("STEP 7: Final Report Generation")
        
        self._generate_final_report(output_dir)
        
        # Step 8: Save configuration and results
        log.section("STEP 8: Save Results")
        
        self._save_results(output_dir)
        
        total_time = time.time() - start_time
        
        log.header("ANALYSIS COMPLETE")
        log.result("Total time", f"{total_time:.1f}", "seconds")
        log.result("Results saved to", output_dir)
        log.result("Mean Binding Energy", f"{mmgbsa_results['mean_binding_energy']:.2f} ± {mmgbsa_results['std_error']:.2f}", "kcal/mol")
        
        return self.results
    
    def _run_entropy_analysis(self, calculator, input_files):
        """
        Run entropy analysis using specified method
        """
        try:
            analysis_settings = self.config['analysis_settings']
            method = analysis_settings.get('entropy_method', 'interaction')
            
            log.info(f"Using entropy method: {method}")
            
            if method == 'interaction':
                # Calculate Interaction Entropy (IE)
                if 'mmgbsa' not in self.results or 'binding_energies' not in self.results['mmgbsa']:
                    log.warning("Interaction Entropy requires binding energies from MM/GBSA run")
                    return None
                
                binding_energies = self.results['mmgbsa']['binding_energies']
                temperature = analysis_settings.get('temperature', 300)
                
                entropy_penalty = calculator.calculate_interaction_entropy(binding_energies, temperature)
                
                log.result("Interaction Entropy (-TΔS)", f"{entropy_penalty:.2f}", "kcal/mol")
                log.info("(Positive value indicates entropic penalty)")
                
                return {
                    'status': 'completed',
                    'method': 'interaction_entropy',
                    'entropy_penalty': entropy_penalty,
                    'T': temperature
                }
                
            elif method in ['nmode', 'normal_mode']:
                # Run Normal Mode Analysis
                log.warning("Running Normal Mode Analysis (Computationally Expensive)")
                
                from .normal_mode import NormalModeAnalysis
                from openmm import unit
                import openmm
                
                try:
                    ligand_mol = self.config['input_files']['ligand_mol']
                    temperature = analysis_settings.get('temperature', 300)
                    
                    log.process("Parameterizing ligand for NMA...")
                    ligand_system, ligand_top, ligand_mol_obj = calculator.parameterize_ligand_openff(ligand_mol)
                    
                    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
                    
                    conf = ligand_mol_obj.conformers[0]
                    if hasattr(conf, 'magnitude'):
                         positions = conf.magnitude * unit.angstrom
                    elif hasattr(conf, 'value_in_unit'):
                         positions = conf.value_in_unit(unit.angstrom) * unit.angstrom
                    else:
                         positions = conf * unit.angstrom
                    
                    log.process("Initializing NormalModeAnalysis...")
                    nma = NormalModeAnalysis(ligand_top, ligand_system, integrator, positions)
                    
                    log.process("Running Minimization (CUDA/CPU)...")
                    nma.CUDAMinimizationCycle(MaxMiniCycle=100)
                    
                    log.process("Calculating Normal Modes...")
                    cutoff_freq = analysis_settings.get('nmode_cutoff_frequency', 10.0)
                    nma.CalculateNormalModes(cutoff_frequency=cutoff_freq)
                    
                    log.process("Calculating Vibrational Entropy...")
                    nma.getVibrationalEntropyQM(Temperature=temperature*unit.kelvin)
                    
                    entropy_val = nma.VibrationalEntropyQM.value_in_unit(unit.kilocalories_per_mole)
                    
                    log.result("Vibrational Entropy Term (TS)", f"{entropy_val:.2f}", "kcal/mol")
                    log.result("Entropic Contribution (-TS)", f"{-entropy_val:.2f}", "kcal/mol")
                    
                    return {
                        'status': 'completed',
                        'method': 'normal_mode',
                        'entropy_term': entropy_val,
                        'contribution': -entropy_val
                    }
                except Exception as e:
                    log.error(f"NMA Failed: {e}")

        except Exception as e:
            log.error(f"Entropy analysis failed: {e}")
            return None
    
    def _run_per_residue_decomposition(self, calculator, input_files, analysis_settings, frame_params=None, output_dir=None):
        """Run per-residue energy decomposition"""
        try:
            decomp_analyzer = PerResidueDecomposition(
                calculator, 
                temperature=analysis_settings.get('temperature', 300),
                output_dir=output_dir
            )
            
            # Extract frame parameters
            if frame_params is None:
                frame_params = {
                    'frame_start': None,
                    'frame_end': None,
                    'frame_stride': None,
                    'frame_selection': 'sequential',
                    'random_seed': 42
                }
            
            decomp_results = decomp_analyzer.run_per_residue_analysis(
                ligand_mol=input_files['ligand_mol'],
                complex_pdb=input_files['complex_pdb'],
                xtc_file=input_files['trajectory'],
                ligand_pdb=input_files['ligand_pdb'],
                max_frames=analysis_settings.get('max_frames', 50),
                decomp_frames=analysis_settings.get('decomp_frames', 10),
                frame_start=frame_params.get('frame_start'),
                frame_end=frame_params.get('frame_end'),
                frame_stride=frame_params.get('frame_stride'),
                frame_selection=frame_params.get('frame_selection', 'sequential'),
                random_seed=frame_params.get('random_seed', 42),
                output_dir=output_dir
            )
            
            return decomp_results
            
        except Exception as e:
            log.error(f"Per-residue decomposition failed: {e}")
            return None
    
    def _generate_final_report(self, output_dir):
        """Generate comprehensive final report"""
        try:
            report_file = output_dir / "final_report.txt"
            
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("MM/GBSA ANALYSIS FINAL REPORT\n")
                f.write("=" * 50 + "\n\n")
                
                # Configuration summary
                f.write("CONFIGURATION SUMMARY:\n")
                f.write("-" * 30 + "\n")
                f.write(f"GB Model: {self.config['analysis_settings'].get('gb_model', 'OBC2')}\n")
                f.write(f"Temperature: {self.config['analysis_settings'].get('temperature', 300)} K\n")
                f.write(f"Salt Concentration: {self.config['analysis_settings'].get('salt_concentration', 0.15)} M\n")
                f.write(f"Max Frames: {self.config['analysis_settings'].get('max_frames', 50)}\n\n")
                
                # MM/GBSA results
                if 'mmgbsa' in self.results:
                    mmgbsa = self.results['mmgbsa']
                    f.write("MM/GBSA RESULTS:\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"Mean Binding Energy: {mmgbsa['mean_binding_energy']:.2f} ± {mmgbsa['std_error']:.2f} kcal/mol\n")
                    f.write(f"Standard Deviation: {mmgbsa['std_dev']:.2f} kcal/mol\n")
                    f.write(f"Frames Analyzed: {mmgbsa['n_frames']}\n")
                    f.write(f"GB Model: {mmgbsa['gb_model']}\n\n")
                    
                    if mmgbsa.get('parameterized_residues'):
                        f.write("PARAMETERIZATION DETAILS:\n")
                        f.write("-" * 25 + "\n")
                        f.write("The following non-standard residues were parameterized on-the-fly via OpenFF:\n")
                        for res in mmgbsa['parameterized_residues']:
                                f.write(f"  • {res}\n")
                        f.write("\n")
                
                # Additional analysis results
                if 'entropy' in self.results:
                    f.write("ENTROPY ANALYSIS:\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"Status: {self.results['entropy']['status']}\n")
                    f.write(f"Method: {self.results['entropy']['method']}\n\n")
                
                if 'decomposition' in self.results:
                    f.write("PER-RESIDUE DECOMPOSITION:\n")
                    f.write("-" * 30 + "\n")
                    f.write("Analysis completed successfully\n\n")
            
            log.success(f"Final report saved: {report_file}")
            
        except Exception as e:
            log.error(f"Report generation failed: {e}")
    
    def _save_results(self, output_dir):
        """Save all results and configuration to output directory"""
        try:
            # Save configuration
            config_file = output_dir / "analysis_config.yaml"
            with open(config_file, 'w', encoding='utf-8') as f:
                yaml.dump(self.config, f, default_flow_style=False, indent=2)
            
            # Save results summary
            results_file = output_dir / "results_summary.yaml"
            with open(results_file, 'w', encoding='utf-8') as f:
                yaml.dump(self.results, f, default_flow_style=False, indent=2)
            
            log.success(f"Configuration saved: {config_file}")
            log.success(f"Results summary saved: {results_file}")
            
        except Exception as e:
            log.error(f"Results saving failed: {e}")


def create_sample_config():
    """Create a sample YAML configuration file"""
    sample_config = {
        'input_files': {
            'ligand_mol': 'test/ligand.sdf',
            'complex_pdb': 'test/complex.pdb',
            'ligand_pdb': 'test/ligand.pdb',
            'trajectory': 'test/complex.xtc'
        },
        'analysis_settings': {
            'temperature': 300,
            'gb_model': 'OBC2',
            'salt_concentration': 0.15,
            'max_frames': 50,
            'verbose': 1,
            'use_cache': True,
            'parallel_processing': False,
            'energy_decomposition': False,
            'run_entropy_analysis': False,
            'run_per_residue_decomposition': False,
            'decomp_frames': 10
        },
        'output_settings': {
            'output_directory': 'mmgbsa_results',
            'save_plots': True,
            'save_trajectories': False
        },
        'advanced_settings': {
            'minimization_tolerance': 1e-6,
            'nma_quality_threshold': 'Good',
            'hot_spot_threshold': -1.0
        }
    }
    
    return sample_config


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='MM/GBSA Analysis Runner')
    parser.add_argument('config_file', nargs='?', help='YAML configuration file')
    parser.add_argument('--create-config', action='store_true', help='Create sample configuration file')
    parser.add_argument('--config-name', default='mmgbsa_config.yaml', help='Name for sample config file')
    
    args = parser.parse_args()
    
    if args.create_config:
        # Create sample configuration
        sample_config = create_sample_config()
        with open(args.config_name, 'w', encoding='utf-8') as f:
            yaml.dump(sample_config, f, default_flow_style=False, indent=2)
        log.success(f"Sample configuration created: {args.config_name}")
        print("Edit this file with your specific parameters and run:")
        print(f"python mmgbsa_runner.py {args.config_name}")
        return
    
    if not args.config_file:
        log.error("Please provide a configuration file or use --create-config")
        print("Usage:")
        print("  python mmgbsa_runner.py config.yaml")
        print("  python mmgbsa_runner.py --create-config")
        return
    
    # Run analysis
    runner = MMGBSARunner(args.config_file)
    results = runner.run_analysis()
    
    if results:
        pass # Already logged by runner
    else:
        log.error("Analysis failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
 