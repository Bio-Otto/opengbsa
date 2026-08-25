"""
Analysis package - modular MM/GBSA analysis pipeline components.

Provides composable analysis engines:
- FrameSelector: Multi-method trajectory frame selection (sequential, equidistant, random)
- EnergyCalculator: Orchestrates binding energy computation across frames
- DecompositionEngine: Per-residue and per-component energy decomposition
- ResultsAggregator: Statistics, validation, convergence analysis
- AnalysisEngine: High-level coordinator delegating to legacy calculator
"""

from typing import Optional, List, Dict, Tuple, Any
import pandas as pd
import numpy as np
from pathlib import Path
import logging

log = logging.getLogger(__name__)


class FrameSelector:
    """
    Multi-method trajectory frame selection for efficient sampling.
    
    Methods:
    - sequential: First N frames (with optional stride)
    - equidistant: Evenly distributed frames across trajectory
    - random: Randomly sampled frames
    """
    
    def __init__(self, trajectory_length: int, verbose: bool = False):
        self.trajectory_length = trajectory_length
        self.verbose = verbose
    
    def select_frames(self, max_frames: Optional[int] = None, frame_start: Optional[int] = None,
                     frame_end: Optional[int] = None, frame_stride: Optional[int] = None,
                     method: str = 'sequential', random_seed: int = 42) -> List[int]:
        """
        Select frames based on parameters and method.
        
        Parameters:
        -----------
        max_frames : int, optional
            Maximum number of frames to select
        frame_start : int, optional
            Start frame (0-indexed)
        frame_end : int, optional
            End frame (0-indexed)
        frame_stride : int, optional
            Frame stride (every Nth frame)
        method : str
            Selection method: 'sequential', 'equidistant', 'random'
        random_seed : int
            Seed for random selection
            
        Returns:
        --------
        list : Selected frame indices
        """
        # Validate frame range
        start = max(0, frame_start if frame_start is not None else 0)
        end = min(self.trajectory_length, frame_end if frame_end is not None else self.trajectory_length)
        
        if start >= end:
            raise ValueError(f"Invalid frame range: {start} >= {end}")
        
        if method == 'sequential':
            frames = self._sequential_selection(start, end, frame_stride, max_frames)
        elif method == 'equidistant':
            frames = self._equidistant_selection(start, end, max_frames)
        elif method == 'random':
            frames = self._random_selection(start, end, frame_stride, max_frames, random_seed)
        else:
            raise ValueError(f"Unknown selection method: {method}")
        
        if self.verbose:
            log.info(f"Selected {len(frames)} frames using '{method}' method from range [{start}, {end})")
        
        return sorted(list(set(frames)))
    
    def _sequential_selection(self, start: int, end: int, stride: Optional[int], 
                             max_frames: Optional[int]) -> List[int]:
        """Sequential frame selection with optional stride."""
        stride = stride or 1
        frames = list(range(start, end, stride))
        if max_frames and len(frames) > max_frames:
            frames = frames[:max_frames]
        return frames
    
    def _equidistant_selection(self, start: int, end: int, max_frames: Optional[int]) -> List[int]:
        """Equidistant frame selection (evenly distributed)."""
        n_frames = end - start
        if max_frames is None or max_frames <= 0:
            max_frames = n_frames
        
        step = n_frames / max(1, max_frames)
        frames = [start + int(i * step) for i in range(max_frames)]
        return [min(f, end - 1) for f in frames]
    
    def _random_selection(self, start: int, end: int, stride: Optional[int], 
                         max_frames: Optional[int], seed: int) -> List[int]:
        """Random frame selection."""
        import random
        random.seed(seed)
        
        available = list(range(start, end, stride or 1))
        if max_frames is None or max_frames <= 0:
            max_frames = len(available)
        
        return random.sample(available, min(max_frames, len(available)))


class EnergyCalculator:
    """
    Orchestrates frame-by-frame energy calculations.
    
    Coordinates:
    - Context preparation (ligand, protein, complex)
    - Per-frame energy extraction
    - Component-wise decomposition (VDW, Electrostatic, GB, SA)
    """
    
    def __init__(self, calculator, verbose: bool = False):
        self.calculator = calculator
        self.verbose = verbose
        self.frame_energies = []
        self.component_energies = {}
    
    def prepare_contexts(self, ligand_system, protein_system, complex_system,
                        ligand_topology, protein_topology, complex_topology,
                        platform, platform_properties: Dict) -> Tuple:
        """
        Prepare OpenMM contexts for energy calculations.
        
        Returns:
        --------
        tuple : (ligand_context, protein_context, complex_context)
        """
        try:
            from openmm import Integrator, VerletIntegrator
            
            # Create dummy integrator (not used, just needed for Context)
            integrator = VerletIntegrator(0.001 * 10**-12)  # dummy timestep
            
            # Create contexts
            ligand_context = platform.createContext(ligand_system, integrator, platform_properties) if platform else None
            protein_context = platform.createContext(protein_system, integrator, platform_properties) if platform else None
            complex_context = platform.createContext(complex_system, integrator, platform_properties) if platform else None
            
            if self.verbose:
                log.info(f"Prepared {sum([1 for c in [ligand_context, protein_context, complex_context] if c])}/3 contexts")
            
            return ligand_context, protein_context, complex_context
        
        except Exception as e:
            if self.verbose:
                log.error(f"Failed to prepare contexts: {e}")
            return None, None, None
    
    def calculate_frame_energy(self, frame_idx: int, positions_dict: Dict,
                              contexts_dict: Dict, indices_dict: Dict) -> Dict[str, float]:
        """
        Calculate binding energy components for a single frame.
        
        Parameters:
        -----------
        frame_idx : int
            Frame index
        positions_dict : dict
            Position arrays: {'ligand': pos, 'protein': pos, 'complex': pos}
        contexts_dict : dict
            OpenMM contexts: {'ligand': ctx, 'protein': ctx, 'complex': ctx}
        indices_dict : dict
            Atom indices: {'ligand': idx, 'protein': idx}
            
        Returns:
        --------
        dict : Energy components {component_name: energy_kcal_mol}
        """
        try:
            from openmm import unit
            
            energy_data = {'frame': frame_idx}
            
            # Set positions and calculate energies
            ligand_ctx = contexts_dict.get('ligand')
            protein_ctx = contexts_dict.get('protein')
            complex_ctx = contexts_dict.get('complex')
            
            ligand_pos = positions_dict.get('ligand')
            protein_pos = positions_dict.get('protein')
            complex_pos = positions_dict.get('complex')
            
            # Calculate component energies
            if complex_ctx and complex_pos is not None:
                complex_ctx.setPositions(complex_pos)
                
                # Total energy
                state = complex_ctx.getState(getEnergy=True)
                energy_data['complex_energy'] = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                
                # Component breakdown (if available)
                # Group 0 (NB): Electrostatics
                try:
                    state_elec = complex_ctx.getState(getEnergy=True, groups={0})
                    energy_data['electrostatic'] = state_elec.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                except:
                    energy_data['electrostatic'] = 0.0
                
                # Group 3 (CustomNB): VDW
                try:
                    state_vdw = complex_ctx.getState(getEnergy=True, groups={3})
                    energy_data['vdw'] = state_vdw.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                except:
                    energy_data['vdw'] = 0.0
                
                # Group 1,2 (GB): Solvation
                try:
                    state_gb = complex_ctx.getState(getEnergy=True, groups={1, 2})
                    energy_data['gb'] = state_gb.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                except:
                    energy_data['gb'] = 0.0
                
                # Group 4 (SA): Surface area
                try:
                    state_sa = complex_ctx.getState(getEnergy=True, groups={4})
                    energy_data['sa'] = state_sa.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                except:
                    energy_data['sa'] = 0.0
            
            # Individual system energies
            if ligand_ctx and ligand_pos is not None:
                ligand_ctx.setPositions(ligand_pos)
                state = ligand_ctx.getState(getEnergy=True)
                energy_data['ligand_energy'] = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            
            if protein_ctx and protein_pos is not None:
                protein_ctx.setPositions(protein_pos)
                state = protein_ctx.getState(getEnergy=True)
                energy_data['protein_energy'] = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            
            # Calculate delta energies
            if 'complex_energy' in energy_data and 'protein_energy' in energy_data and 'ligand_energy' in energy_data:
                energy_data['binding_energy'] = (energy_data['complex_energy'] - 
                                               energy_data['protein_energy'] - 
                                               energy_data['ligand_energy'])
            else:
                energy_data['binding_energy'] = 0.0
            
            return energy_data
        
        except Exception as e:
            if self.verbose:
                log.error(f"Error calculating frame {frame_idx}: {e}")
            return {'frame': frame_idx, 'error': str(e)}
    
    def calculate_frame_energies(self, frames: List[int], trajectory_file: str, 
                                topology_file: str, contexts_dict: Dict,
                                indices_dict: Dict,
                                output_csv: Optional[str] = None) -> pd.DataFrame:
        """
        Calculate binding energies for selected frames.
        
        Parameters:
        -----------
        frames : list
            Frame indices to analyze
        trajectory_file : str
            Path to trajectory file
        topology_file : str
            Path to topology file
        contexts_dict : dict
            OpenMM contexts for calculations
        indices_dict : dict
            Atom selection indices
        output_csv : str, optional
            Output CSV file for results
            
        Returns:
        --------
        pd.DataFrame : Energy results for each frame
        """
        try:
            import mdtraj as md
            
            results = []
            traj = md.load(trajectory_file, top=topology_file)
            
            for i, frame_idx in enumerate(frames):
                if self.verbose and (i + 1) % 10 == 0:
                    log.info(f"Processed {i + 1}/{len(frames)} frames")
                
                # Extract frame
                if frame_idx >= len(traj):
                    if self.verbose:
                        log.warning(f"Frame {frame_idx} out of range (trajectory length: {len(traj)})")
                    continue
                
                frame = traj[frame_idx]
                
                # Build position dictionaries (placeholder - actual values from legacy)
                positions_dict = {
                    'complex': frame.xyz[0] * 10,  # Convert nm to angstrom for OpenMM
                    'ligand': None,  # Would be extracted from indices
                    'protein': None
                }
                
                # Calculate energy
                energy_data = self.calculate_frame_energy(frame_idx, positions_dict, 
                                                         contexts_dict, indices_dict)
                results.append(energy_data)
                self.frame_energies.append(energy_data)
            
            df = pd.DataFrame(results)
            if output_csv:
                df.to_csv(output_csv, index=False)
                if self.verbose:
                    log.info(f"Energy results saved to {output_csv}")
            
            return df
        
        except Exception as e:
            if self.verbose:
                log.error(f"Failed to calculate frame energies: {e}")
            return pd.DataFrame()


class DecompositionEngine:
    """
    Handles per-residue and per-component energy decomposition.
    
    Provides:
    - Per-residue energy contributions
    - Interaction analysis (residue-ligand)
    - Hot spot identification
    """
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.decomposition_data = []
        self.per_residue_energies = {}
    
    def decompose_per_residue(self, results_df: pd.DataFrame, 
                             decomp_data: Optional[Dict] = None,
                             decomp_type: str = 'full') -> pd.DataFrame:
        """
        Decompose binding energy into per-residue contributions.
        
        Parameters:
        -----------
        results_df : pd.DataFrame
            Frame-level energy results
        decomp_data : dict, optional
            Pre-computed decomposition data (from legacy calculator)
        decomp_type : str
            Decomposition type: 'full' (all contributions), 'mmgbsa' (MM-PBSA std),
            'simplified' (hot-spots only)
            
        Returns:
        --------
        pd.DataFrame : Per-residue decomposition results
        """
        try:
            if decomp_data is None:
                # If no data provided, create template structure
                decomp_data = {}
            
            per_residue_list = []
            
            if decomp_type == 'full':
                # Full decomposition includes all energy terms
                per_residue_list = self._full_decomposition(results_df, decomp_data)
            
            elif decomp_type == 'mmgbsa':
                # Standard MM-PBSA style decomposition
                per_residue_list = self._mmgbsa_decomposition(results_df, decomp_data)
            
            elif decomp_type == 'simplified':
                # Hot-spot residues only (top contributors)
                per_residue_list = self._simplified_decomposition(results_df, decomp_data)
            
            df_decomp = pd.DataFrame(per_residue_list)
            
            if self.verbose and len(df_decomp) > 0:
                log.info(f"Decomposed into {len(df_decomp)} residue entries")
                if 'total' in df_decomp.columns:
                    top_contribs = df_decomp.nsmallest(5, 'total')
                    log.info("Top 5 most favorable residues:")
                    for _, row in top_contribs.iterrows():
                        log.info(f"  {row.get('residue_name', '?')}{row.get('residue_number', '?')}: {row.get('total', 0):.2f} kcal/mol")
            
            self.decomposition_data = per_residue_list
            return df_decomp
        
        except Exception as e:
            if self.verbose:
                log.error(f"Decomposition failed: {e}")
            return pd.DataFrame()
    
    def _full_decomposition(self, results_df: pd.DataFrame, decomp_data: Dict) -> List[Dict]:
        """Full energy decomposition per residue."""
        results = []
        
        # Template structure
        if 'residues' in decomp_data:
            for res_info in decomp_data['residues']:
                result_entry = {
                    'residue_name': res_info.get('name', 'UNK'),
                    'residue_number': res_info.get('number', 0),
                    'chain': res_info.get('chain', 'X'),
                    'vdw': res_info.get('vdw', 0.0),
                    'electrostatic': res_info.get('elec', 0.0),
                    'gb': res_info.get('gb', 0.0),
                    'sa': res_info.get('sa', 0.0),
                }
                result_entry['total'] = sum([
                    result_entry.get('vdw', 0),
                    result_entry.get('electrostatic', 0),
                    result_entry.get('gb', 0),
                    result_entry.get('sa', 0)
                ])
                results.append(result_entry)
        
        return results
    
    def _mmgbsa_decomposition(self, results_df: pd.DataFrame, decomp_data: Dict) -> List[Dict]:
        """MM-PBSA style decomposition."""
        # Simplified MM-PBSA: VDW + Elec + GB + SA
        return self._full_decomposition(results_df, decomp_data)
    
    def _simplified_decomposition(self, results_df: pd.DataFrame, decomp_data: Dict) -> List[Dict]:
        """Simplified hot-spot decomposition (top N contributors)."""
        full_decomp = self._full_decomposition(results_df, decomp_data)
        
        # Sort by total contribution and take top N
        sorted_decomp = sorted(full_decomp, key=lambda x: x.get('total', 0))[:10]
        
        if self.verbose:
            log.info(f"Identified {len(sorted_decomp)} hot-spot residues")
        
        return sorted_decomp
    
    def identify_hotspot_residues(self, decomp_df: pd.DataFrame, threshold: float = -1.0) -> List[str]:
        """
        Identify hot-spot residues with significant contributions.
        
        Parameters:
        -----------
        decomp_df : pd.DataFrame
            Decomposition results
        threshold : float
            Energy threshold for hot-spot (kcal/mol)
            
        Returns:
        --------
        list : Residue identifiers of hot-spots
        """
        if 'total' not in decomp_df.columns:
            return []
        
        hotspots = decomp_df[decomp_df['total'] < threshold].copy()
        
        result = []
        for _, row in hotspots.iterrows():
            res_id = f"{row.get('residue_name', '?')}{row.get('residue_number', '?')}"
            result.append(res_id)
        
        return result


class ResultsAggregator:
    """Aggregates, validates, and analyzes results."""
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.statistics = {}
    
    def compute_statistics(self, binding_energies: np.ndarray) -> Dict[str, float]:
        """Compute comprehensive binding energy statistics."""
        stats = {
            'mean': np.mean(binding_energies),
            'std_dev': np.std(binding_energies),
            'median': np.median(binding_energies),
            'q1': np.percentile(binding_energies, 25),
            'q3': np.percentile(binding_energies, 75),
            'min': np.min(binding_energies),
            'max': np.max(binding_energies),
            'std_error': np.std(binding_energies) / np.sqrt(len(binding_energies))
        }
        return stats
    
    def bootstrap_uncertainty(self, binding_energies: np.ndarray, 
                            n_bootstrap: int = 1000) -> Dict[str, float]:
        """Calculate uncertainty using bootstrap resampling."""
        np.random.seed(42)
        bootstrap_means = []
        n = len(binding_energies)
        
        for _ in range(n_bootstrap):
            sample = np.random.choice(binding_energies, size=n, replace=True)
            bootstrap_means.append(np.mean(sample))
        
        bootstrap_means = np.array(bootstrap_means)
        return {
            'mean': np.mean(bootstrap_means),
            'std': np.std(bootstrap_means),
            'ci_lower': np.percentile(bootstrap_means, 2.5),
            'ci_upper': np.percentile(bootstrap_means, 97.5)
        }
    
    def validate_results(self, results_df: pd.DataFrame) -> List[str]:
        """Validate MM/GBSA results for reasonableness."""
        warnings = []
        
        if 'binding_energy' not in results_df.columns:
            return warnings
        
        be = results_df['binding_energy']
        mean_be = be.mean()
        std_be = be.std()
        
        if mean_be > 10:
            warnings.append(f"Binding energy very positive ({mean_be:.1f} kcal/mol) - check for errors")
        elif mean_be < -50:
            warnings.append(f"Binding energy very negative ({mean_be:.1f} kcal/mol) - check for errors")
        
        if std_be > 10:
            warnings.append(f"High standard deviation ({std_be:.1f} kcal/mol) - system may be unstable")
        
        if len(be) >= 20:
            first_half = be[:len(be)//2].mean()
            second_half = be[len(be)//2:].mean()
            if abs(first_half - second_half) > 2.0:
                warnings.append(f"Poor convergence: halves differ by {abs(first_half - second_half):.1f} kcal/mol")
        
        return warnings
    
    def check_convergence(self, binding_energies: np.ndarray, 
                         window_size: int = 10) -> Dict[str, Any]:
        """Check if binding energy has converged."""
        if len(binding_energies) < 2 * window_size:
            return {'converged': False, 'reason': 'insufficient_data'}
        
        df = pd.Series(binding_energies)
        running_avg = df.rolling(window=window_size).mean()
        recent_avg = running_avg.tail(window_size).mean()
        early_avg = running_avg.iloc[window_size:2*window_size].mean()
        
        threshold = 1.0  # kcal/mol
        converged = abs(recent_avg - early_avg) < threshold
        
        return {
            'converged': converged,
            'threshold': threshold,
            'difference': abs(recent_avg - early_avg)
        }


class AnalysisEngine:
    """
    High-level analysis coordinator -- currently a thin pass-through.

    `run`/`run_comprehensive` both delegate directly to the underlying
    `calculator` (in practice `mmgbsa.mmgbsa_core.GBSACalculator`, the actual
    engine that performs frame selection, energy calculation, decomposition,
    and result aggregation internally).

    NOTE: this class also constructs `self.frame_selector`, `self.energy_calc`,
    `self.decomposition`, and `self.aggregator` (the modular-refactor
    components below), but none of them are actually called by `run` or
    `run_comprehensive` -- they are inert scaffolding from an in-progress
    refactor, not yet wired into the execution path. Do not assume using
    `AnalysisEngine` gets you the modular pipeline; all real computation
    currently happens inside `calculator`.
    """

    def __init__(self, calculator=None, verbose: bool = False):
        self.verbose = verbose
        if calculator is None:
            from .. import mmgbsa_core
            self.calculator = mmgbsa_core.GBSACalculator(verbose=verbose)
        else:
            self.calculator = calculator

        self.frame_selector = FrameSelector(100, verbose=verbose)  # Will be updated with actual length
        self.energy_calc = EnergyCalculator(self.calculator, verbose=verbose)
        self.decomposition = DecompositionEngine(verbose=verbose)
        self.aggregator = ResultsAggregator(verbose=verbose)

    def run(self, *args, **kwargs):
        """Run analysis by delegating directly to `self.calculator.run`
        (the sibling components constructed above are not used here)."""
        return self.calculator.run(*args, **kwargs)

    def run_comprehensive(self, *args, **kwargs):
        """
        Run the comprehensive MM/GBSA pipeline with validation and analysis.
        
        Coordinates:
        1. Frame selection (sequential/equidistant/random)
        2. Energy calculations
        3. Results validation
        4. Bootstrap uncertainty
        5. Convergence analysis
        6. Report generation
        """
        if hasattr(self.calculator, 'run_comprehensive'):
            return self.calculator.run_comprehensive(*args, **kwargs)
        else:
            return self.calculator.run(*args, **kwargs)
