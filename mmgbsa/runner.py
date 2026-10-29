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
from .mmgbsa_core import GBSACalculator
from .interaction_entropy import describe, entropy_result
from .decomposition import PerResidueDecomposition
from .logger import ToolLogger
from .reporting import HTMLReportGenerator
from .config import ConfigManager, OPTIONS, effective_config

# Initialize global logger
log = ToolLogger()

def make_run_directory(base_dir, prefix='analysis'):
    """Create and return a new folder <base_dir>/<prefix>_<YYYYmmdd_HHMMSS>.

    A run started in the same second as another gets the suffix _2, _3, ...
    The folder is created with mkdir without exist_ok, so two processes
    (e.g. a job array) never share one.
    """
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    n = 1
    while True:
        path = base / (f"{prefix}_{stamp}" if n == 1 else f"{prefix}_{stamp}_{n}")
        try:
            path.mkdir()
            return path
        except FileExistsError:
            n += 1


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
            Path to a YAML configuration file, validated as the opengbsa
            command validates it, or a configuration dict (not validated).
        output_dir : str, optional
            Folder of this run, used as given. Without it the run creates a new
            folder <output_directory>/analysis_<date>_<time>/.
        """
        if isinstance(config_file, (str, Path)):
            self.config_file = str(config_file)
            config = self._load_config()
        else:
            self.config_file = None
            config = config_file
        # Every option present, with its default where the configuration leaves it
        # out (mmgbsa.config.OPTIONS); written to analysis_config.yaml
        self.config = effective_config(config)
        self.results = {}
        self.output_dir = output_dir

    def _load_config(self):
        """Load and validate a YAML configuration file (as the opengbsa command does)."""
        manager = ConfigManager()
        if not manager.load_config(self.config_file):
            raise ValueError(f"Could not read the configuration file {self.config_file}")
        if not manager.validate_config():
            raise ValueError("Configuration validation failed:\n  - "
                             + "\n  - ".join(manager.get_validation_errors()))
        for warning in manager.get_validation_warnings():
            log.warning(warning)
        return manager.get_config()

    def _validate_input_files(self):
        """Every input file the configuration names must exist."""
        input_files = self.config['input_files']
        missing_files = []
        if not input_files.get('complex_pdb'):
            missing_files.append("complex_pdb: not given")
        for key, option in OPTIONS['input_files'].items():
            file_path = input_files.get(key)
            if not option.path or file_path is None:
                continue
            if not Path(str(file_path)).exists():
                missing_files.append(f"{key}: {file_path}")

        if missing_files:
            log.error("Missing input files:")
            for missing in missing_files:
                log.error(f"  • {missing}")
            return False

        log.success("All required input files found")
        return True

    def _create_output_directory(self):
        """The run folder: output_dir as given, else a new timestamped folder."""
        if self.output_dir is None:
            base = self.config.get('output_settings', {}).get('output_directory') or 'mmgbsa_results'
            output_path = make_run_directory(base, 'analysis')
        else:
            output_path = Path(self.output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

        # Update instance variable so it can be retrieved
        self.output_dir = output_path

        log.info(f"Output directory established: {output_path}")
        return output_path

    def _setup_mmgbsa_calculator(self):
        """Setup MM/GBSA calculator with configuration"""
        analysis_settings = self.config['analysis_settings']
        forcefield_settings = self.config.get('forcefield_settings', {})

        calculator = GBSACalculator(
            temperature=analysis_settings.get('temperature', 300),
            verbose=bool(analysis_settings.get('verbose', True)),
            gb_model=analysis_settings.get('gb_model', 'OBC2'),
            salt_concentration=analysis_settings.get('salt_concentration', 0.15),
            max_workers=analysis_settings.get('max_workers', None),
            protein_forcefield=forcefield_settings.get('protein_forcefield', 'amber'),
            charge_method=analysis_settings.get('charge_method', 'am1bcc'),
            solute_dielectric=analysis_settings.get('solute_dielectric', 1.0),
            solvent_dielectric=analysis_settings.get('solvent_dielectric', None),   # 78.5, or MMPBSA.py's with amber_compatible
            entropy_method=(analysis_settings.get('entropy_method', 'interaction')
                            if analysis_settings.get('run_entropy_analysis', False) else 'none'),
            decomposition_method=analysis_settings.get('decomposition_method', 'full'),
            sa_model=analysis_settings.get('sa_model', None),                       # LCPO
            reimage_trajectory=analysis_settings.get('reimage_trajectory', True),
            charmm_params=forcefield_settings.get('charmm_params', None),
            charmm_coordinates=forcefield_settings.get('charmm_coordinates', None),
            gb_radii_set=analysis_settings.get('gb_radii_set', 'auto'),
            amber_compatible=analysis_settings.get('amber_compatible', False),
            rgbmax=analysis_settings.get('rgbmax', None),
            lcpo_typing=analysis_settings.get('lcpo_typing', None),     # 'element'; 'amber' with amber_compatible
            nmode_frames=analysis_settings.get('nmode_frames', 10),       # null: every analysed frame
            nmode_drms=analysis_settings.get('nmode_drms', 0.001),
            nmode_maxcyc=analysis_settings.get('nmode_maxcyc', 10000),
            nmode_precision=analysis_settings.get('nmode_precision', 'mixed'),
        )

        # Apply platform separation settings (analysis/decomposition) when provided.
        platform_settings = self.config.get('platform_settings', {})
        if platform_settings and hasattr(calculator, 'set_platform_settings'):
            calculator.set_platform_settings(platform_settings)
        
        return calculator
    
    def _resolve_binding_mode_selections(self):
        """
        Translate `binding_mode` config option into receptor/ligand MDTraj selections.
        
        Modes:
          standard -- receptor = everything but the ligand, ligand detected (or
            ligand_selection). (The deprecated 'dimer_ligand' is mapped to it by
            mmgbsa.config.)
          ppi -- chain 0 = receptor, chain 1 = ligand; the complex must then have
            exactly two chains (checked by run_comprehensive).
        Explicit receptor_selection/ligand_selection always win. With only
        receptor_selection the ligand is every other atom.
        
        Returns:
            (receptor_selection, ligand_selection, require_two_chains); the
            selections may be None (= auto-detect).
        """
        analysis_settings = self.config.get('analysis_settings', {})
        binding_mode = analysis_settings.get('binding_mode', 'standard')
        
        # Manual overrides always win
        manual_rec = analysis_settings.get('receptor_selection')
        manual_lig = analysis_settings.get('ligand_selection')
        require_two_chains = False
        
        if binding_mode == 'ppi':
            receptor_selection = manual_rec if manual_rec else 'chainid 0'
            ligand_selection   = manual_lig if manual_lig else 'chainid 1'
            require_two_chains = not (manual_rec and manual_lig)
            log.info(f"Binding mode: ppi  (receptor='{receptor_selection}', ligand='{ligand_selection}')")
        else:
            receptor_selection = manual_rec  # None = auto (all non-ligand)
            ligand_selection   = manual_lig  # None = auto (detect by resname)
            if receptor_selection and not ligand_selection:
                ligand_selection = f"not ({receptor_selection})"
                receptor_selection = None
                log.info(f"Binding mode: standard  (ligand = every atom not in receptor_selection "
                         f"'{manual_rec}')")
            elif ligand_selection and not receptor_selection:
                log.info(f"Binding mode: standard  (ligand='{ligand_selection}', receptor = every other atom)")
            elif ligand_selection:
                log.info(f"Binding mode: standard  (receptor='{receptor_selection}', ligand='{ligand_selection}')")
            else:
                log.info("Binding mode: standard  (receptor = everything but the ligand, ligand detected)")
        
        return receptor_selection, ligand_selection, require_two_chains

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
        self.calculator = calculator
        
        # Step 4: Run MM/GBSA analysis
        log.section("STEP 4: MM/GBSA Analysis")
        
        input_files = self.config['input_files']
        analysis_settings = self.config['analysis_settings']
        
        # A .pdb complex is Coordinate Mode. GROMACS Native Mode is chosen by giving
        # the .top as complex_pdb; a .top next to a .pdb no longer switches the mode.
        complex_pdb = input_files.get('complex_pdb')

        # GROMACS AUTO-CONVERSION CHECK
        if complex_pdb and str(complex_pdb).endswith('.top'):
            log.info("Detected GROMACS topology (.top). Attempting auto-conversion to Amber format...")
            try:
                from .conversion import GromacsPreprocessor
                from pathlib import Path
                
                top_path = Path(complex_pdb)
                # Infer coordinate file (.gro)
                gro_files = list(top_path.parent.glob("*.gro"))
                coord_file = None
                
                # Logic to find best matching .gro
                if gro_files:
                     for g in gro_files:
                         if g.stem == top_path.stem:
                             coord_file = g
                             break
                     if not coord_file:
                         coord_file = gro_files[0]
                
                if coord_file:
                    log.info(f"Found GROMACS coordinate file: {coord_file}")
                    conversion_results = GromacsPreprocessor.convert_to_amber(
                        top_path, coord_file, output_dir=Path(self.output_dir) / "gromacs_conversion",
                        gb_model=analysis_settings.get('gb_model', 'OBC2')
                    )
                    
                    # UPDATE INPUT FILES WITH CONVERTED TOPOLOGIES
                    log.success("GROMACS conversion successful. Updating input configuration.")
                    input_files['complex_pdb'] = conversion_results['topology']       # dry complex
                    input_files['receptor_topology'] = conversion_results['receptor_topology']
                    input_files['ligand_topology'] = conversion_results['ligand_topology']
                    # We MUST use the .gro file as the solvated topology for MDTraj. 
                    # If we use the generated .prmtop, we lose the experimental PDB residue numbers 
                    # because PRMTOP format natively drops resSeq strings and forces a 1..N order.
                    input_files['solvated_topology'] = str(coord_file)
                    
                else:
                    log.warning("Could not find .gro file for GROMACS conversion. Proceeding without conversion (may fail).")
            except Exception as e:
                log.error(f"GROMACS auto-conversion failed: {e}")
                # We continue, likely to fail later, but error log is imperative.
        
        # OPENMM SYSTEM (.xml + .pdb of the same name): convert to an Amber topology,
        # which then runs in Native Mode like any .prmtop
        if complex_pdb and str(complex_pdb).endswith('.xml'):
            from pathlib import Path   # Path is local to this method (imported in branches above)
            from .conversion import OpenMMSystemPreprocessor
            pdb_file = OpenMMSystemPreprocessor.find_pdb(complex_pdb)
            log.info(f"Detected a serialised OpenMM System; converting it with {pdb_file.name} "
                     f"to an Amber topology...")
            input_files['complex_pdb'] = OpenMMSystemPreprocessor.convert_to_amber(
                complex_pdb, pdb_file, Path(self.output_dir) / "openmm_system",
                gb_model=analysis_settings.get('gb_model', 'OBC2'),
                radii_set=analysis_settings.get('gb_radii_set', 'auto'))
            complex_pdb = input_files['complex_pdb']
            if not input_files.get('solvated_topology'):
                # the PDB keeps the residue numbering of the system for the trajectory
                input_files['solvated_topology'] = str(pdb_file)
            log.success(f"OpenMM System converted: {input_files['complex_pdb']}")

        # Get frame selection parameters for reporting
        frame_params = {
            'max_frames': analysis_settings.get('max_frames', 50),
            'frame_start': analysis_settings.get('frame_start', None),
            'frame_end': analysis_settings.get('frame_end', None),
            'frame_stride': analysis_settings.get('frame_stride', None),
            'frame_selection': analysis_settings.get('frame_selection', 'sequential')
        }
        
        # Resolve receptor/ligand selections from binding_mode (or manual overrides)
        resolved_receptor_sel, resolved_ligand_sel, require_two_chains = self._resolve_binding_mode_selections()
        
        mmgbsa_results = calculator.run_comprehensive(
            ligand_mol=input_files.get('ligand_mol'),
            complex_pdb=input_files['complex_pdb'],
            xtc_file=input_files.get('trajectory'),
            ligand_pdb=input_files.get('ligand_pdb'),
            max_frames=analysis_settings.get('max_frames', 50),
            energy_decomposition=analysis_settings.get('energy_decomposition', False),
            frame_start=frame_params.get('frame_start'),
            frame_end=frame_params.get('frame_end'),
            frame_stride=frame_params.get('frame_stride'),
            frame_selection=frame_params.get('frame_selection', 'sequential'),
            random_seed=analysis_settings.get('random_seed', 42),
            qha_analyze_complex=analysis_settings.get('qha_analyze_complex', False),
            output_dir=output_dir,
            ligand_selection=resolved_ligand_sel,
            receptor_selection=resolved_receptor_sel,
            receptor_topology=input_files.get('receptor_topology'),
            ligand_topology=input_files.get('ligand_topology'),
            solvated_topology=input_files.get('solvated_topology'),
            print_interval=analysis_settings.get('print_interval', 10),
            require_two_chains=require_two_chains,
        )
        
        if not mmgbsa_results:
            log.error("MM/GBSA analysis failed!")
            return None
        
        self.results['mmgbsa'] = mmgbsa_results
        log.success(f"MM/GBSA (Standard) completed: {mmgbsa_results['mean_binding_energy']:.2f} ± {mmgbsa_results['std_error']:.2f} kcal/mol")

        if mmgbsa_results.get('parameterized_residues'):
            log.result("Parameterized Residues", ", ".join(mmgbsa_results['parameterized_residues']))

        # Step 5: Run entropy analysis (if enabled)
        if analysis_settings.get('run_entropy_analysis', False):
            log.section("STEP 5: Entropy Analysis")
            
            entropy_results = self._run_entropy_analysis(calculator, input_files)
            if entropy_results:
                self.results['entropy'] = entropy_results
                log.success("Entropy analysis completed")
        
        # Step 6: Run per-residue decomposition (if enabled)
        decomposition_failure = None
        if analysis_settings.get('run_per_residue_decomposition', False):
            log.section("STEP 6: Per-Residue Decomposition")
            
            # Get frame parameters for per-residue decomposition
            frame_params = {
                'frame_start': analysis_settings.get('frame_start', None),
                'frame_end': analysis_settings.get('frame_end', None),
                'frame_stride': analysis_settings.get('frame_stride', None),
                'frame_selection': analysis_settings.get('frame_selection', 'sequential'),
                'random_seed': analysis_settings.get('random_seed', 42)
            }
            
            decomp_results, decomp_error = self._run_per_residue_decomposition(
                calculator, input_files, analysis_settings, frame_params, output_dir,
                receptor_selection=resolved_receptor_sel,
                ligand_selection=resolved_ligand_sel,
                mmgbsa_results=mmgbsa_results
            )
            if decomp_error is None:
                decomp_results['status'] = 'completed'
                self.results['decomposition'] = decomp_results
                log.success("Per-residue decomposition completed")
            else:
                # never reported as success: marked in results_summary.yaml and
                # final_report.txt; with strict_platform the run stops (below)
                self.results['decomposition'] = {'status': 'failed', 'error': decomp_error}
                decomposition_failure = decomp_error
                log.error(f"PER-RESIDUE DECOMPOSITION FAILED: {decomp_error}. The MM/GBSA results "
                          f"are complete; the per-residue results were not written.")
        
        # Step 7: Generate final report
        log.section("STEP 7: Final Report Generation")
        
        self._generate_final_report(output_dir)
        
        # Step 8: Save configuration and results
        log.section("STEP 8: Save Results")
        
        self._save_results(output_dir)

        if decomposition_failure and (self.config.get('platform_settings') or {}).get('strict_platform'):
            raise RuntimeError(f"Per-residue decomposition failed ({decomposition_failure}); "
                               f"platform_settings.strict_platform: true stops the run with an error. "
                               f"The MM/GBSA results are in {output_dir}.")

        total_time = time.time() - start_time
        
        log.header("ANALYSIS COMPLETE")
        log.result("Total time", f"{total_time:.1f}", "seconds")
        log.result("Results saved to", str(output_dir))
        
        # Report Standard Deviation (SD) instead of SEM
        mean_val = mmgbsa_results.get('mean_binding_energy', 0.0)
        std_dev = mmgbsa_results.get('std_dev', 0.0)
        if 'std_dev' not in mmgbsa_results and 'binding_energies' in mmgbsa_results:
             import numpy as np
             std_dev = np.std(mmgbsa_results['binding_energies'])
             
        log.result("Mean Binding Energy", f"{mean_val:.2f} ± {std_dev:.2f}", "kcal/mol (SD)")

        # Sanity-check the final result against known physically-reasonable
        # ranges (implausibly large |dG|, unstable/high-variance runs) so
        # problems surface immediately rather than silently making it into a
        # results table. `components`-based checks (vdw/gb sign sanity) are
        # not wired here: `mmgbsa_results` does not expose mean component
        # values under a stable key, only per-frame arrays on the calculator
        # instance (`self.energies['complex_gb']` etc.), so passing a guessed
        # mapping would silently never fire rather than actually validate.
        try:
            from .validation import TopologyValidator
            sanity_warnings = TopologyValidator.validate_system_sanity(
                binding_energy=mean_val, std_dev=std_dev, components=None
            )
            for w in sanity_warnings:
                if w.severity == 'ERROR':
                    log.error(str(w))
                else:
                    log.warning(str(w))
        except Exception as e:
            log.warning(f"Result sanity-check skipped due to an internal error: {e}")

        return self.results
    
    def _run_entropy_analysis(self, calculator, input_files):
        """
        Run entropy analysis using specified method
        """
        try:
            analysis_settings = self.config['analysis_settings']
            method = analysis_settings.get('entropy_method', 'interaction')
            
            log.info(f"Using entropy method: {method}")
            
            if method in ('interaction', 'c2'):
                # Interaction entropy (IE) or its second-order cumulant (C2),
                # both from the gas-phase receptor-ligand interaction energy
                if 'mmgbsa' not in self.results or 'interaction_energies' not in self.results['mmgbsa']:
                    log.warning("Entropy estimates require the per-frame interaction energies from the MM/GBSA run")
                    return None

                temperature = analysis_settings.get('temperature', 300)
                entropy = entropy_result(self.results['mmgbsa']['interaction_energies'], method, temperature)

                for line in describe(entropy):
                    log.info(line)
                log.result(f"{entropy['label']} (-TΔS)", f"{entropy['value']:.2f}", "kcal/mol")
                log.info("(Positive value indicates entropic penalty)")

                mmgbsa = self.results['mmgbsa']
                if not mmgbsa.get('entropy'):
                    # Report it like an entropy_method result (HTML report, dG)
                    mmgbsa['entropy'] = entropy
                    mmgbsa['entropy_penalty'] = entropy['value']
                    mmgbsa['delta_g'] = mmgbsa['mean_binding_energy'] + entropy['value']

                return {
                    'status': 'completed',
                    'method': entropy['label'],
                    'entropy_penalty': entropy['value'],
                    'T': temperature,
                    **{k: v for k, v in entropy.items() if k not in ('method', 'label', 'value')}
                }
                
            elif method == 'quasiharmonic':
                # Quasi-harmonic entropy, computed by the calculator from the
                # covariance of the analysed frames (mmgbsa.quasi_harmonic)
                entropy = self.results.get('mmgbsa', {}).get('entropy')
                if not entropy or entropy.get('method') != 'quasiharmonic':
                    log.warning("Quasi-harmonic entropy is not available from the MM/GBSA run")
                    return None
                if entropy.get('value') is None:
                    log.result("Ligand TS (absolute, quasi-harmonic)", f"{entropy['TS_ligand']:.2f}", "kcal/mol")
                    log.info("No binding entropy: set qha_analyze_complex: true for -TdS of binding")
                else:
                    log.result(f"{entropy['label']} (-TΔS)", f"{entropy['value']:.2f}", "kcal/mol")
                return {
                    'status': 'completed',
                    'method': entropy['label'],
                    'entropy_penalty': entropy.get('value'),
                    'T': entropy['temperature'],
                    **{k: v for k, v in entropy.items() if k not in ('method', 'label', 'value')}
                }

            elif method in ['nmode', 'normal_mode']:
                # Normal-mode binding entropy, computed by the calculator from the
                # minimized complex, receptor and ligand (mmgbsa.nmode_entropy)
                entropy = self.results.get('mmgbsa', {}).get('entropy')
                if not entropy or entropy.get('method') != 'nmode':
                    log.warning("Normal-mode entropy is not available from the MM/GBSA run")
                    return None
                log.result(f"{entropy['label']} (-TΔS)", f"{entropy['value']:.2f} ± {entropy['sd']:.2f}",
                           "kcal/mol")
                return {
                    'status': 'completed',
                    'method': entropy['label'],
                    'entropy_penalty': entropy['value'],
                    'T': entropy['temperature'],
                    **{k: v for k, v in entropy.items() if k not in ('method', 'label', 'value')}
                }

        except Exception as e:
            log.error(f"Entropy analysis failed: {e}")
            return None
    
    def _run_per_residue_decomposition(self, calculator, input_files, analysis_settings, frame_params=None, output_dir=None,
                                        receptor_selection=None, ligand_selection=None, mmgbsa_results=None):
        """Run per-residue energy decomposition"""
        try:
            # Enable Parallel Processing
            # Get n_jobs from advanced settings, default to -1 (All Cores) for maximum performance
            n_jobs = self.config.get('advanced_settings', {}).get('n_jobs', -1)
            
            decomp_analyzer = PerResidueDecomposition(
                calculator, 
                temperature=analysis_settings.get('temperature', 300),
                output_dir=output_dir,
                n_jobs=n_jobs,
                report_raw_energies=analysis_settings.get('report_raw_energies', False)
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
            
            # When custom selections triggered distillation in run_comprehensive, the
            # calculator stores the distilled paths as instance attributes.  Use those
            # so decomposition analyses exactly the same atom subset as the main calculation.
            decomp_complex_pdb = getattr(calculator, 'distilled_complex_pdb', None) or input_files['complex_pdb']
            decomp_xtc_file    = getattr(calculator, 'distilled_xtc_file',    None) or input_files['trajectory']
            if hasattr(calculator, 'distilled_complex_pdb'):
                log.info(f"Decomposition: using distilled complex PDB: {decomp_complex_pdb}")
            # Files reduced to the selected atoms (distilled.xtc) hold no solvent, so the
            # solvated topology no longer matches them; a Native-Mode complex reduced to
            # the selected atoms takes its parameters from the reduced files, since the
            # input topology's atom indices no longer apply
            distilled = bool(getattr(calculator, 'distilled_files', None))
            decomp_solvated = None if distilled else input_files.get('solvated_topology')
            decomp_parameters = (decomp_complex_pdb if getattr(calculator, 'native_reduced', False)
                                 else input_files.get('complex_pdb'))

            # The decomposed frames: a subset of the analysed frames (mmgbsa.frames)
            decomposition_frames = None
            try:
                decomposition_frames = calculator.decomposition_frames(
                    analysis_settings.get('decomp_frames', 10),
                    frame_params.get('frame_selection', 'sequential'), frame_params.get('random_seed', 42))
                log.info(f"Decomposition frames ({len(decomposition_frames)}): "
                         f"{decomposition_frames[:10]}{' ...' if len(decomposition_frames) > 10 else ''}")
            except Exception as e:
                log.warning(f"Could not determine the decomposition frames: {e}")

            decomp_results = decomp_analyzer.run_per_residue_analysis(
                ligand_mol=input_files.get('ligand_mol'),
                complex_pdb=decomp_complex_pdb,
                xtc_file=decomp_xtc_file,
                ligand_pdb=input_files.get('ligand_pdb'),
                max_frames=analysis_settings.get('max_frames', 50),
                decomp_frames=analysis_settings.get('decomp_frames', 10),
                frame_start=frame_params.get('frame_start'),
                frame_end=frame_params.get('frame_end'),
                frame_stride=frame_params.get('frame_stride'),
                frame_selection=frame_params.get('frame_selection', 'sequential'),
                random_seed=frame_params.get('random_seed', 42),
                solvated_topology=decomp_solvated,
                receptor_topology=input_files.get('receptor_topology'),
                ligand_topology=input_files.get('ligand_topology'),
                ligand_resname=getattr(calculator, 'ligand_resname', input_files.get('ligand_resname')),
                decomposition_topology_file=decomp_parameters,
                output_dir=output_dir,
                plot_top_residues=analysis_settings.get('plot_top_residues', 10),
                receptor_selection=getattr(calculator, 'receptor_selection', receptor_selection),
                ligand_selection=getattr(calculator, 'ligand_selection', ligand_selection),
                skip_baseline=True,
                baseline_results=mmgbsa_results
            )

            
            if not decomp_results or decomp_results.get('decomposition_results') is None:
                # run_per_residue_analysis returns the MM/GBSA results when the
                # decomposition itself failed
                return None, ("the decomposition produced no per-residue results "
                              "(see the errors logged above)")
            decomp_results['frames'] = decomposition_frames

            # Generate Interactive HTML Report
            if decomp_results and hasattr(decomp_analyzer, 'frame_data'):
                try:
                    # Select best PDB for visualization
                    # Prioritize the temporary PDB extracted from trajectory if available
                    viz_pdb = input_files.get('complex_pdb')
                    temp_pdb = output_dir / "temp_from_traj.pdb"
                    if temp_pdb.exists():
                        viz_pdb = str(temp_pdb)
                    
                    pandamap_path = output_dir / "structure_3d.html"
                    # IMPORTANT: HTML reporter expects an HTML file for iframe embedding.
                    # Do not fall back to a PDB path here, otherwise raw PDB text is embedded
                    # into srcdoc and breaks interactive report rendering.
                    if pandamap_path.exists():
                        pandamap_path = str(pandamap_path.resolve())
                    else:
                        pandamap_path = None

                    html_gen = HTMLReportGenerator(output_dir, config=self.config)
                    html_gen.generate_report(
                        analysis_results=decomp_results,
                        frame_data=decomp_analyzer.frame_data,
                        global_results=decomp_results.get('mmgbsa_results'),
                        complex_pdb_path=pandamap_path,
                        ligand_resname=input_files.get('ligand_resname') or self.config.get('input_files', {}).get('ligand_resname'),
                        decomposition_frames=decomposition_frames
                    )
                except Exception as e:
                    log.warning(f"Interactive report generation failed: {e}")
            
            return decomp_results, None

        except Exception as e:
            import traceback
            log.error(f"Per-residue decomposition failed: {e}")
            traceback.print_exc()
            return None, f"{type(e).__name__}: {e}"
    
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
                        
                    f.write("PHYSICS ASSUMPTIONS AND FALLBACKS:\n")
                    f.write("-" * 35 + "\n")
                    if not mmgbsa.get('physics_assumptions'):
                        f.write("  • None. The pipeline utilized explicitly structured parameters and pure explicit configurations without defaulting to empirical heuristics.\n\n")
                    else:
                        seen = set()
                        for assumption in mmgbsa['physics_assumptions']:
                            if assumption not in seen:
                                f.write(f"  • {assumption}\n")
                                seen.add(assumption)
                        f.write("\n")
                
                # Additional analysis results
                if 'entropy' in self.results:
                    f.write("ENTROPY ANALYSIS:\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"Status: {self.results['entropy']['status']}\n")
                    f.write(f"Method: {self.results['entropy']['method']}\n")
                    if 'sigma_interaction_energy' in self.results['entropy']:
                        for line in describe(self.results['entropy']):
                            f.write(f"{line}\n")
                    if 'TS_ligand' in self.results['entropy']:
                        e = self.results['entropy']
                        f.write(f"Ligand TS = {e['TS_ligand']:.2f} kcal/mol (quasi-harmonic, "
                                f"{e.get('n_frames', 0)} frames at {e['T']:.2f} K)\n")
                        if e.get('entropy_penalty') is not None:
                            f.write(f"-TdS (binding) = {e['entropy_penalty']:.2f} kcal/mol\n")
                        else:
                            f.write("No binding entropy (qha_analyze_complex: false)\n")
                    if 'n_converged' in self.results['entropy']:
                        e = self.results['entropy']
                        f.write(f"-TdS = {e['entropy_penalty']:.2f} +/- {e['sd']:.2f} kcal/mol (SD) over "
                                f"{e['n_converged']} of {e['n_frames']} frames at {e['T']:.2f} K\n")
                        f.write(f"  translational {-e['TdS_translational']:.2f}, rotational "
                                f"{-e['TdS_rotational']:.2f}, vibrational {-e['TdS_vibrational']:.2f} kcal/mol\n")
                    f.write("\n")
                
                if 'decomposition' in self.results:
                    decomposition = self.results['decomposition']
                    f.write("PER-RESIDUE DECOMPOSITION:\n")
                    f.write("-" * 30 + "\n")
                    if decomposition.get('status') == 'failed':
                        f.write(f"FAILED: {decomposition.get('error')}\n\n")
                    else:
                        frames = decomposition.get('frames')
                        f.write("Completed" + (f" on {len(frames)} frames" if frames else "") + "\n\n")
            
            log.success(f"Final report saved: {report_file}")
            
        except Exception as e:
            log.error(f"Report generation failed: {e}")
    
    def _save_results(self, output_dir):
        """Save all results and configuration to output directory"""
        try:
            # Save the effective configuration: every option with the value this
            # run used, so the run can be repeated from the file
            config_file = output_dir / "analysis_config.yaml"
            with open(config_file, 'w', encoding='utf-8') as f:
                f.write("# Effective configuration of this run: every option with the value used\n"
                        "# (defaults filled in). Repeat the run with: opengbsa analysis_config.yaml\n")
                yaml.safe_dump(self._sanitize_results_for_yaml(self.effective_run_config()), f,
                               default_flow_style=False, indent=2, sort_keys=False)
            
            # Save results summary
            results_file = output_dir / "results_summary.yaml"
            with open(results_file, 'w', encoding='utf-8') as f:
                yaml.safe_dump(self._sanitize_results_for_yaml(self.results), f, default_flow_style=False, indent=2, sort_keys=False)
            
            log.success(f"Configuration saved: {config_file}")
            log.success(f"Results summary saved: {results_file}")
            
        except Exception as e:
            log.error(f"Results saving failed: {e}")

    def effective_run_config(self):
        """The configuration with the values the run resolved at run time.

        Options whose null means "decided by the run" get the value used:
        max_workers, sa_model, lcpo_typing, rgbmax and solvent_dielectric (which
        amber_compatible takes from the topology). Input paths are those used,
        after any format conversion.
        """
        import copy
        config = copy.deepcopy(self.config)
        calculator = getattr(self, 'calculator', None)
        if calculator is not None:
            settings = config['analysis_settings']
            settings['max_workers'] = calculator.max_workers
            settings['sa_model'] = calculator.sa_model
            settings['lcpo_typing'] = calculator.lcpo_typing
            settings['rgbmax'] = calculator.rgbmax
            if calculator.solvent_dielectric is not None:
                settings['solvent_dielectric'] = float(calculator.solvent_dielectric)
        return config

    def _sanitize_results_for_yaml(self, data):
        """
        Convert runtime result objects (DataFrame, ndarray, Path, numpy scalars, etc.)
        into compact YAML-safe summaries.
        """
        try:
            import pandas as pd
            import numpy as np
        except Exception:
            pd = None
            np = None

        if data is None:
            return None

        # Convert numpy scalar types before primitive checks; np.float64 is also
        # an instance of float and would otherwise leak into safe_dump unchanged.
        if np is not None and isinstance(data, np.generic):
            return data.item()
        if isinstance(data, (str, int, float, bool)):
            return data
        if isinstance(data, Path):
            return str(data)
        if np is not None and isinstance(data, np.ndarray):
            return {
                'type': 'ndarray',
                'shape': list(data.shape),
                'dtype': str(data.dtype)
            }
        if pd is not None and isinstance(data, pd.DataFrame):
            return {
                'type': 'DataFrame',
                'rows': int(len(data)),
                'columns': [str(c) for c in data.columns.tolist()]
            }

        if isinstance(data, dict):
            return {str(k): self._sanitize_results_for_yaml(v) for k, v in data.items()}
        if isinstance(data, (list, tuple, set)):
            out = [self._sanitize_results_for_yaml(v) for v in list(data)]
            if len(out) > 200:
                return {
                    'type': type(data).__name__,
                    'count': len(out),
                    'sample': out[:50]
                }
            return out

        return str(data)


def main():
    """`python -m mmgbsa.runner config.yaml`: the same run as `opengbsa config.yaml`."""
    parser = argparse.ArgumentParser(description='MM/GBSA Analysis Runner')
    parser.add_argument('config_file', nargs='?', help='YAML configuration file')
    parser.add_argument('--create-config', action='store_true', help='Create sample configuration file')
    parser.add_argument('--config-name', default='mmgbsa_config.yaml', help='Name for sample config file')

    args = parser.parse_args()

    if args.create_config:
        ConfigManager().create_default_config(args.config_name)
        log.success(f"Sample configuration created: {args.config_name}")
        print(f"Edit this file with your specific parameters and run: opengbsa {args.config_name}")
        return

    if not args.config_file:
        log.error("Please provide a configuration file or use --create-config")
        print("Usage:")
        print("  python -m mmgbsa.runner config.yaml")
        print("  python -m mmgbsa.runner --create-config")
        sys.exit(1)

    try:
        runner = MMGBSARunner(args.config_file)
    except ValueError as e:
        log.error(str(e))
        sys.exit(1)
    results = runner.run_analysis()
    if not results:
        log.error("Analysis failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
