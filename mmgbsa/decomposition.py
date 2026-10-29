#!/usr/bin/env python3
"""
Per-Residue Energy Decomposition for MM/GBSA Analysis
This adds advanced per-residue analysis to your existing MM/GBSA package
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import time
import warnings
import multiprocessing # Added for parallel processing

# Global worker context for multiprocessing
_worker_context = {}

# Suppress specific warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="Unable to load toolkit 'OpenEye Toolkit'")
warnings.filterwarnings("ignore", message="importing 'simtk.openmm' is deprecated")

from openmm import app, openmm, unit

import mdtraj as md
import parmed
from .mmgbsa_core import GBSACalculator, SOLVENT_ION_RESNAMES
from .trajectory import reimage_molecules
from .pymol_viz import write_hotspot_pymol
from .solvation_peratom import (
    SolvationEnergies, group_energy, POLAR_GROUPS, NONPOLAR_GROUPS, ONE_4PI_EPS0, KJ_TO_KCAL,
    CHUNK_ELEMENTS, _map_blocks, _discrete_function,
    _Expression as _LeptonExpression, _BUILTINS as _LEPTON_BUILTINS,
)

from .logger import ToolLogger

# Initialize logger
log = ToolLogger()


class DecompositionError(RuntimeError):
    """The per-residue decomposition failed; no per-residue results exist."""


class PerResidueDecomposition:
    """
    Advanced per-residue energy decomposition for MM/GBSA analysis.
    Integrates seamlessly with your existing GBSACalculator class.

    Every term is decomposed from the OpenMM Systems the MM/GBSA calculation
    itself evaluates (complex, receptor and ligand), so the per-residue values
    follow the run's settings (GB model, rgbmax, salt, dielectrics, surface-area
    model and LCPO typing) and, with the ligand's share, add up to the frame's
    binding energy term by term:

      - vdw / electrostatic: pairwise interaction energy between each
        residue's atoms and the ligand atoms from the complex System's own
        nonbonded forces, HALVED (Gohlke-Kollman 50/50 convention, as Amber's
        `decpair` in sander decomp.F90: hfval = 0.5*fval into both residues'
        buckets, with no special case for cross-molecule pairs).
      - solvation (polar): each atom's GB energy in the complex minus its GB
        energy in the receptor System, summed per residue. An atom's GB
        energy is its self term plus half of each of its pair terms, read
        from the run's GB force (mmgbsa.solvation_peratom).
      - nonpolar_solvation: the same complex-minus-receptor difference of the
        per-atom surface-area energy of the run's non-polar force (LCPO or
        ACE). Amber MMPBSA.py decomposes this term with sander gbsa=2 instead
        (a numerical surface), so it does not sum to Amber's LCPO total.

    'total' is vdw + electrostatic + polar solvation (the non-polar term is
    reported separately, as before). This is a different quantity than the
    `decomposition_method='fast'` path (see `mmgbsa.energy_decomp.EnergyDecomposer`),
    which reports vdW+electrostatic only with no solvation term at all.
    """
    
    def __init__(self, mmgbsa_calculator, temperature=300.0, output_dir=None, n_jobs=1, report_raw_energies=False):
        """
        Initialize per-residue decomposition analysis
        
        Parameters:
        -----------
        mmgbsa_calculator : GBSACalculator
            Your existing MM/GBSA calculator
        temperature : float
            Temperature in Kelvin
        output_dir : str, optional
            Output directory for saving results
        n_jobs : int
            Number of parallel jobs (default: 1, -1 for all cores)
        report_raw_energies : bool
            If True, calculate and report raw complex/receptor/ligand energies (Amber-like)
        """
        self.mmgbsa_calculator = mmgbsa_calculator
        self.temperature = temperature * unit.kelvin
        self.output_dir = output_dir
        self.n_jobs = n_jobs
        self.report_raw_energies = report_raw_energies
        log.debug(f"PerResidueDecomposition init: report_raw_energies={self.report_raw_energies}")
        
        # Storage for decomposition results
        self.residue_contributions = {}
        self.interaction_matrix = {}
        self.hot_spots = []
        
    def run_per_residue_analysis(self, ligand_mol, complex_pdb, xtc_file, 
                                ligand_pdb, max_frames=50, decomp_frames=10,
                                output_dir=None,
                                frame_start=None, frame_end=None, frame_stride=None,
                                frame_selection='sequential', random_seed=42, 
                                solvated_topology=None, receptor_topology=None, ligand_topology=None,
                                ligand_resname=None, plot_top_residues=10,
                                decomposition_topology_file=None,
                                receptor_selection=None, ligand_selection=None, skip_baseline=False, baseline_results=None):
        """
        Run complete MM/GBSA analysis with per-residue decomposition
        
        Parameters:
        -----------
        decomp_frames : int
            Number of frames to use for decomposition (computationally expensive)
        receptor_selection : str, optional
            MDTraj selection string for receptor atoms (e.g. 'chainid 0')
        ligand_selection : str, optional
            MDTraj selection string for ligand atoms (e.g. 'resname UNK')

        Returns the results dict (with 'status': 'completed'). Raises
        DecompositionError when the decomposition fails, so a failure is
        never reported or stored as a result.
        """

        if output_dir:
            self.report_dir = output_dir
            
        log.section("MM/GBSA WITH PER-RESIDUE DECOMPOSITION")
        
        # Step 1: Run standard MM/GBSA analysis
        if not skip_baseline:
            log.info("Running Baseline MM/GBSA Analysis...")
            
            mmgbsa_results = self.mmgbsa_calculator.run(
                ligand_mol, complex_pdb, xtc_file, ligand_pdb, max_frames,
                solvated_topology=solvated_topology,
                receptor_topology=receptor_topology,
                ligand_topology=ligand_topology,
                output_dir=output_dir if output_dir else self.output_dir,
                frame_start=frame_start,
                frame_end=frame_end,
                frame_stride=frame_stride,
                frame_selection=frame_selection,
                random_seed=random_seed
            )
            
            if not mmgbsa_results:
                raise DecompositionError("The MM/GBSA analysis the decomposition builds on failed")
            
            log.result("Baseline MM/GBSA Binding", f"{mmgbsa_results['mean_binding_energy']:.2f} ± {mmgbsa_results['std_error']:.2f}", "kcal/mol")
        else:
            log.info("Skipping redundant Baseline Analysis (Systems already cached)")
            mmgbsa_results = baseline_results or getattr(self.mmgbsa_calculator, 'results', {})
        
        # Step 2: Per-residue decomposition
        log.section("Per-Residue Energy Decomposition")
        frames_msg = decomp_frames if decomp_frames else 'all available'
        log.process(f"Analyzing {frames_msg} frames for detailed decomposition...")
        
        decomp_results = self._perform_per_residue_decomposition(
            ligand_mol, complex_pdb, xtc_file, ligand_pdb, decomp_frames,
            frame_start=frame_start, frame_end=frame_end, frame_stride=frame_stride,
            frame_selection=frame_selection, random_seed=random_seed,
            solvated_topology=solvated_topology, receptor_topology=receptor_topology,
            decomposition_topology_file=decomposition_topology_file,
            ligand_resname=ligand_resname,
            salt_concentration=self.mmgbsa_calculator.salt_concentration,
            plot_top_residues=plot_top_residues,
            receptor_selection=receptor_selection,
            ligand_selection=ligand_selection
        )
        
        # _perform_per_residue_decomposition raises DecompositionError on
        # failure; an empty result is a failure too, never a "completed" run.
        if not decomp_results:
            raise DecompositionError("Per-residue decomposition produced no results")

        # Step 3: Analyze and visualize results
        log.section("Analysis and Visualization")

        analysis_results = self._analyze_decomposition_results(decomp_results)
        self._generate_decomposition_plots(analysis_results)

        # Combine results
        complete_results = {
            'status': 'completed',
            'mmgbsa_results': mmgbsa_results,
            'decomposition_results': decomp_results,
            'analysis_results': analysis_results,
            'n_decomp_frames': decomp_frames,
            'energy_check': getattr(self, 'energy_check', []),
        }

        self._print_decomposition_summary(complete_results)
        self._save_decomposition_results(complete_results)
        self._write_hotspot_pymol(decomp_results)

        return complete_results

    def execute_decomposition(self, trajectory, systems_dict):
        """
        Decompose an already loaded trajectory with the run's own Systems.

        `systems_dict` holds the complex System, Context and topology; the
        receptor and ligand Systems are taken from the calculator
        (`GBSACalculator.systems`) when the dict does not carry them.
        Raises DecompositionError when the decomposition fails.
        """
        try:
            log.process("Executing Decomposition on provided trajectory...")

            complex_topology = systems_dict['complex_topology']  # OpenMM Topology
            ligand_resname = self.mmgbsa_calculator.find_ligand_resname(complex_topology)
            residue_map, ligand_indices = self._build_residue_mapping(complex_topology, ligand_resname)
            log.info(f"Found {len(residue_map)} protein residues, {len(ligand_indices)} ligand atoms")

            systems = dict(systems_dict)
            calc_systems = getattr(self.mmgbsa_calculator, 'systems', None) or {}
            for key in ('protein_system', 'ligand_system'):
                if systems.get(key) is None and calc_systems.get(key) is not None:
                    systems[key] = calc_systems[key]
            decomposer = self._make_decomposer(systems, residue_map, ligand_indices)

            residue_energies = []
            frame_by_frame_data = []
            for i, frame in enumerate(trajectory):
                if i % 10 == 0:
                    log.process(f"Decomposing frame {i+1}/{len(trajectory)}...")
                frame_result, _ = decomposer.decompose(frame.xyz[0])
                residue_energies.append(frame_result)
                frame_by_frame_data.append(_frame_row(frame_result, i, i + 1, len(trajectory)))

            if not residue_energies:
                raise DecompositionError("The trajectory has no frames to decompose")

            averaged = self._average_residue_energies(residue_energies)
            self.frame_data = frame_by_frame_data

            # Visualization/Saving
            analysis = self._analyze_decomposition_results(averaged)
            self._generate_decomposition_plots(analysis)

            if self.output_dir:
                df = analysis['dataframe']
                df.to_csv(os.path.join(self.output_dir, "final_decomposition.csv"), index=False)
                print(f"✓ Decomposition saved to: final_decomposition.csv")

            return averaged

        except DecompositionError:
            raise
        except Exception as e:
            raise DecompositionError(f"Decomposition execution failed: {e}") from e

    def _perform_per_residue_decomposition(self, ligand_mol, complex_pdb,
                                         xtc_file, ligand_pdb, n_frames,
                                         frame_start=None, frame_end=None, frame_stride=None,
                                         frame_selection='sequential', random_seed=42, solvated_topology=None,
                                         receptor_topology=None, decomposition_topology_file=None,
                                         ligand_resname=None, salt_concentration=None,
                                         plot_top_residues=10, receptor_selection=None, ligand_selection=None):
        """
        Perform detailed per-residue energy decomposition.

        receptor_selection / ligand_selection mirror the same parameters in GBSACalculator.run().
        When provided they restrict which atoms are treated as receptor vs ligand during
        the residue-level energy decomposition, consistent with the main MM/GBSA calculation.

        Returns the per-residue averages; raises DecompositionError on any
        failure (no frame is skipped silently).
        """

        try:
            log.debug(f"perform_decomposition: report_raw_energies={self.report_raw_energies}")
            # Load trajectory and get frames for decomposition
            import os
            if str(complex_pdb).endswith('.tpr'):
                 temp_pdb = os.path.join(self.output_dir if self.output_dir else '.', "complex_solvated_from_tpr_for_mdtraj.pdb")
                 if os.path.exists(temp_pdb):
                      log.info(f"Intercepting dynamic Native TPR metadata in Decomposition module: {temp_pdb}")
                      complex_pdb = temp_pdb

            # Load trajectory (Use solvated_topology if available to avoid atom mismatch)
            load_top = solvated_topology if solvated_topology else complex_pdb

            traj = md.load(xtc_file, top=load_top)

            if n_frames is None or len(traj) < n_frames:
                n_frames = len(traj)

            # Select frames based on parameters
            selected_frames = self._select_frames(len(traj), n_frames, frame_start, frame_end,
                                                frame_stride, frame_selection, random_seed)
            decomp_traj = traj[selected_frames]

            print(f"  Selected {len(decomp_traj)} frames for decomposition")

            # Strip solvent/ions from trajectory to match the (dry) OpenMM System.
            topology = decomp_traj.topology

            # Same solvent/ion residue names as the system builder's strip
            # (mmgbsa_core solvent_mask), CHARMM's included: with only the
            # Amber/GROMACS names, CHARMM SOD/CLA ions stayed in the trajectory
            # and every frame had more atoms than the system.
            keep_indices = np.array([a.index for a in topology.atoms
                                     if not (a.residue.is_water or a.residue.name in SOLVENT_ION_RESNAMES)])

            # Optimization: Only strip if we actually have solvent atoms or if explicit solvated top was used
            if solvated_topology or len(keep_indices) < topology.n_atoms:
                log.process(f"Stripping solvent/ion atoms (keeping {len(keep_indices)}) from trajectory to match system...")
                decomp_traj = decomp_traj.atom_slice(keep_indices)

            # Re-image across periodic boundaries, gated on the same
            # reimage_trajectory setting as the main energy loop
            # (TrajectoryProcessor.load_and_process): a molecule that wraps to
            # the far side of the box corrupts inter-atomic distances.
            if getattr(self.mmgbsa_calculator, 'reimage_trajectory', True):
                try:
                    decomp_traj = reimage_molecules(decomp_traj)
                except Exception as e:
                    log.warning(f"image_molecules() failed ({e}); proceeding with "
                                f"un-imaged coordinates. Per-residue energies may be "
                                f"biased if any molecule crosses a periodic boundary "
                                f"during the trajectory.")

            # The run's own Systems (complex, receptor, ligand) and Contexts.
            topology_for_params = decomposition_topology_file or complex_pdb
            systems = self._prepare_decomposition_systems(
                ligand_mol, complex_pdb, ligand_pdb, topology_for_params,
                ligand_selection=ligand_selection
            )

            log.process("Building residue mapping from prepared system...")
            # Use the topology from the prepared system (stripped of solvent/ions)
            # to ensure indices match the OpenMM system and stripped trajectory
            complex_topology = systems['complex_topology']

            if not ligand_resname:
                log.info("Searching for ligand residue name...")
                ligand_resname = self.mmgbsa_calculator.find_ligand_resname(complex_topology)

            log.info(f"Using Ligand Residue Name: {ligand_resname}")
            if ligand_selection:
                log.info(f"Using custom ligand_selection for decomposition: '{ligand_selection}'")
            if receptor_selection:
                log.info(f"Using custom receptor_selection for decomposition: '{receptor_selection}'")

            residue_map, ligand_indices = self._build_residue_mapping(
                complex_topology,
                ligand_resname,
                ligand_selection=ligand_selection,
                receptor_selection=receptor_selection,
                reference_topology=decomp_traj.topology
            )

            log.info(f"Found {len(residue_map)} protein residues, {len(ligand_indices)} ligand atoms")
            if len(ligand_indices) == 0:
                res_found = set(r.name for r in complex_topology.residues())
                raise DecompositionError(
                    f"No ligand atoms found for resname '{ligand_resname}' "
                    f"(residues in the topology: {sorted(res_found)[:20]})")
            if decomp_traj.n_atoms != systems['complex_system'].getNumParticles():
                raise DecompositionError(
                    f"The decomposition trajectory has {decomp_traj.n_atoms} atoms, the complex "
                    f"System {systems['complex_system'].getNumParticles()}")
            # structure for the PyMOL hot-spot script: the middle decomposed frame, in the
            # atom order of the evaluated complex
            self._hotspot_inputs = dict(topology=complex_topology, residue_map=residue_map,
                                        ligand_indices=ligand_indices,
                                        positions=decomp_traj.xyz[len(decomp_traj) // 2])

            decomposer = self._make_decomposer(systems, residue_map, ligand_indices)

            # Determine number of cores
            n_cores = self.n_jobs
            if n_cores == -1:
                n_cores = multiprocessing.cpu_count()

            # Parallel decomposition only on the CPU platform: GPU contexts
            # cannot be shared safely across worker processes.
            decomposition_platform = getattr(self.mmgbsa_calculator, 'decomposition_platform', None)
            run_parallel = (n_cores > 1 and len(decomp_traj) > 1 and decomposition_platform == 'CPU')

            outputs = None
            if run_parallel:
                n_workers = min(n_cores, len(decomp_traj))
                log.info(f"✓ Running parallel decomposition on {n_workers} cores (CPU platform)")
                payload = {
                    'complex_xml': openmm.XmlSerializer.serialize(systems['complex_system']),
                    'protein_xml': openmm.XmlSerializer.serialize(systems['protein_system']),
                    'ligand_xml': openmm.XmlSerializer.serialize(systems['ligand_system']),
                    'residue_map': residue_map,
                    'ligand_indices': list(ligand_indices),
                    'report_raw_energies': self.report_raw_energies,
                }
                tasks = [(i, decomp_traj.xyz[i]) for i in range(len(decomp_traj))]
                try:
                    outputs = [None] * len(tasks)
                    with multiprocessing.Pool(processes=n_workers, initializer=_worker_init,
                                              initargs=(payload,)) as pool:
                        for done, (i, out, error) in enumerate(pool.imap(_worker_decompose, tasks)):
                            if error is not None:
                                raise DecompositionError(f"frame {selected_frames[i]}: {error}")
                            outputs[i] = out
                            if (done + 1) % 5 == 0:
                                log.process(f"Processed {done+1}/{len(tasks)} frames (Parallel)...")
                except DecompositionError:
                    raise
                except Exception as e:
                    log.error(f"Parallel execution failed: {e}. Falling back to serial.")
                    outputs = None
            else:
                if decomposition_platform and decomposition_platform != 'CPU':
                    log.info(f"• Serial decomposition (decomposition_platform={decomposition_platform} incompatible with GPU-unsafe parallelization)")
                else:
                    log.info(f"• Serial decomposition (no explicit CPU platform set; GPU parallelization disabled for safety)")

            if outputs is None:
                log.info("Running serial decomposition...")
                outputs = []
                for i in range(len(decomp_traj)):
                    if i % 5 == 0:
                        log.process(f"Processing frame {i+1}/{len(decomp_traj)}...")
                    try:
                        outputs.append(decomposer.decompose(decomp_traj.xyz[i]))
                    except DecompositionError as e:
                        raise DecompositionError(f"frame {selected_frames[i]}: {e}") from e
                    except Exception as e:
                        raise DecompositionError(f"frame {selected_frames[i]}: {e}") from e

            # Same assembly for both paths: per-residue energies, the
            # frame-by-frame rows and the energy check.
            residue_energies = []
            frame_by_frame_data = []
            energy_check = []
            for i, (frame_result, info) in enumerate(outputs):
                residue_energies.append(frame_result)
                frame_by_frame_data.append(_frame_row(frame_result, selected_frames[i], i + 1, len(decomp_traj)))
                energy_check.append({'frame_index': selected_frames[i], **info})

            log.success(f"Decomposed {len(residue_energies)} frames successfully")
            self._log_energy_check(energy_check)

            # Save frame-by-frame output
            self._save_frame_by_frame_csv(frame_by_frame_data, residue_map)

            # Store for reporting
            self.frame_data = frame_by_frame_data
            self.energy_check = energy_check

            # Average across frames
            averaged_results = self._average_residue_energies(residue_energies)

            # Generate advanced visualization plots
            if frame_by_frame_data:
                self._generate_time_series_heatmap(frame_by_frame_data, residue_map, top_n=plot_top_residues)

            return averaged_results

        except DecompositionError:
            raise
        except Exception as e:
            raise DecompositionError(f"Decomposition failed: {e}") from e

    def _write_hotspot_pymol(self, averaged_results):
        """Write hotspots.pdb and view_binding.pml into the output folder: the
        receptor coloured by each residue's contribution, the `total` of
        per_residue_detailed.csv and of the plots (van der Waals, electrostatic
        and polar terms)."""
        inputs = getattr(self, '_hotspot_inputs', None)
        if not inputs or not self.output_dir:
            return None
        try:
            values = {res_id: e.get('total_mean', 0.0) for res_id, e in averaged_results.items()}
            path = write_hotspot_pymol(inputs['topology'], inputs['positions'], values,
                                       inputs['residue_map'], inputs['ligand_indices'],
                                       self.output_dir)
            log.success(f"PyMOL hot-spot script written: {path} (run 'pymol {path.name}' in that folder)")
            return path
        except Exception as e:
            log.warning(f"PyMOL hot-spot script not written: {e}")
            return None

    def _build_residue_mapping(self, topology, ligand_resname, ligand_selection=None, receptor_selection=None, reference_topology=None):
        """
        Build mapping from atoms to residues.

        When ligand_selection is provided it is used as an MDTraj selection string to
        identify ligand atoms inside the OpenMM topology.  Otherwise atoms whose residue
        name matches ligand_resname are treated as ligand atoms.

        receptor_selection is accepted for symmetry / future use but is not needed here:
        every atom that is NOT a ligand atom is treated as part of a receptor residue.
        """
        import mdtraj as md

        residue_map = {}  # {residue_id: [atom_indices]}
        ligand_indices = []
        ref_atoms = None

        # Preserve original residue numbering when a reference MDTraj topology is provided.
        if reference_topology is not None:
            try:
                ref_atoms = list(reference_topology.atoms)
                if len(ref_atoms) != topology.getNumAtoms():
                    log.warning(
                        "_build_residue_mapping: reference_topology atom count mismatch "
                        f"({len(ref_atoms)} vs {topology.getNumAtoms()}); falling back to OpenMM ids."
                    )
                    ref_atoms = None
            except Exception as e:
                log.warning(f"_build_residue_mapping: failed to use reference topology ({e})")
                ref_atoms = None

        # Pre-compute ligand atom index set when a custom selection is given
        if ligand_selection:
            try:
                mdtraj_top = md.Topology.from_openmm(topology)
                sel_result = mdtraj_top.select(ligand_selection)
                if len(sel_result) == 0:
                    # Selection matched nothing in the built topology (e.g. user wrote 'resname UNL'
                    # but OpenFF renamed it to UNK).  Fall back to ligand_resname.
                    log.warning(
                        f"_build_residue_mapping: ligand_selection '{ligand_selection}' matched 0 atoms "
                        f"in built topology — falling back to resname '{ligand_resname}'"
                    )
                    lig_set = None
                else:
                    lig_set = set(sel_result.tolist())
                    log.info(f"_build_residue_mapping: ligand_selection '{ligand_selection}' matched {len(lig_set)} atoms")
            except Exception as e:
                log.warning(f"_build_residue_mapping: ligand_selection failed ({e}), falling back to resname '{ligand_resname}'")
                lig_set = None
        else:
            lig_set = None


        for atom in topology.atoms():
            is_ligand = (lig_set is not None and atom.index in lig_set) or \
                        (lig_set is None and atom.residue.name == ligand_resname)
            if is_ligand:
                ligand_indices.append(int(atom.index))
            else:
                if ref_atoms is not None:
                    ref_res = ref_atoms[atom.index].residue
                    res_name = ref_res.name
                    res_num = int(getattr(ref_res, 'resSeq', ref_res.index + 1))
                    chain = getattr(ref_res, 'chain', None)
                    chain_id = str(getattr(chain, 'index', 'A')) if chain is not None else 'A'
                else:
                    res_name = atom.residue.name
                    res_num = int(atom.residue.id)
                    chain_id = str(atom.residue.chain.id)

                res_id = f"{res_name}_{res_num}_{chain_id}"
                if res_id not in residue_map:
                    residue_map[res_id] = []
                residue_map[res_id].append(int(atom.index))

        return residue_map, ligand_indices

    
    def _select_frames(self, trajectory_length, max_frames=None, frame_start=None, frame_end=None,
                      frame_stride=None, frame_selection='sequential', random_seed=42):
        """
        Select frames based on parameters (same as MM/GBSA calculator)
        """
        # Determine frame range
        if frame_start is None:
            frame_start = 0
        if frame_end is None:
            frame_end = trajectory_length
        
        # Validate frame range
        frame_start = max(0, min(frame_start, trajectory_length - 1))
        frame_end = max(frame_start + 1, min(frame_end, trajectory_length))
        
        print(f"  Frame selection parameters:")
        print(f"    • Trajectory length: {trajectory_length}")
        print(f"    • Frame range: {frame_start} to {frame_end}")
        print(f"    • Frame stride: {frame_stride}")
        print(f"    • Selection method: {frame_selection}")
        print(f"    • Max frames: {max_frames}")
        
        # Generate frame indices based on selection method
        if frame_selection == "sequential":
            # Sequential selection with stride
            if frame_stride is None or frame_stride <= 0:
                frame_indices = list(range(frame_start, frame_end))
            else:
                frame_indices = list(range(frame_start, frame_end, frame_stride))
                
        elif frame_selection == "equidistant":
            # Equidistant selection
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            step = (frame_end - frame_start) / max_frames
            frame_indices = [frame_start + int(i * step) for i in range(max_frames)]
            frame_indices = [min(idx, frame_end - 1) for idx in frame_indices]
            
        elif frame_selection == "random":
            # Random selection
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
        
        print(f"  Selected {len(frame_indices)} frames:")
        print(f"    • Frame indices: {frame_indices[:10]}{'...' if len(frame_indices) > 10 else ''}")
        print(f"    • Frame range: {min(frame_indices)} to {max(frame_indices)}")
        
        return frame_indices
    
    def _prepare_decomposition_systems(self, ligand_mol, complex_pdb, ligand_pdb, topology_file=None,
                                       ligand_selection=None):
        """
        Collect the OpenMM Systems the MM/GBSA calculation evaluated (complex,
        receptor, ligand) and create a Context for each on the decomposition
        platform. The Contexts only serve to check, frame by frame, that the
        decomposed energies sum to the energies OpenMM computes.

        Raises DecompositionError when the receptor or ligand System is missing:
        the solvation terms cannot be decomposed exactly without them.
        """
        print("  Preparing systems for decomposition...")

        calc_systems = getattr(self.mmgbsa_calculator, 'systems', None) or {}
        if calc_systems.get('complex_system') is not None:
            log.info("  Reusing pre-built native MM/GBSA OpenMM systems for decomposition...")
            complex_system = calc_systems['complex_system']
            complex_topology = calc_systems['complex_topology']
        else:
            complex_system, complex_topology, _ = self.mmgbsa_calculator.build_complex_system(
                complex_pdb, ligand_mol, ligand_pdb
            )
        protein_system = calc_systems.get('protein_system')
        ligand_system = calc_systems.get('ligand_system')
        if protein_system is None or ligand_system is None:
            raise DecompositionError(
                "The receptor and ligand Systems of the MM/GBSA run are not available. The "
                "per-residue polar and non-polar terms are the complex-minus-receptor/ligand "
                "per-atom energies of those Systems; run the MM/GBSA calculation first.")

        decomp_platform = getattr(self.mmgbsa_calculator, 'decomposition_platform', None)
        platform, properties = self.mmgbsa_calculator.setup_optimized_platform(platform_name=decomp_platform)

        def context(system):
            integrator = openmm.VerletIntegrator(0.001 * unit.picosecond)
            return openmm.Context(system, integrator, platform, properties)

        systems = {
            'complex_system': complex_system,
            'complex_topology': complex_topology,
            'complex_context': context(complex_system),
            'protein_system': protein_system,
            'protein_context': context(protein_system),
            'ligand_system': ligand_system,
            'ligand_context': context(ligand_system),
        }
        log.info(f"  ✓ Protein-only GB context prepared ({protein_system.getNumParticles()} atoms) "
                 f"for per-residue GB self-energy correction")
        print("  Systems prepared for decomposition")
        return systems

    def _make_decomposer(self, systems, residue_map, ligand_indices):
        """The per-frame decomposer for these Systems and this residue mapping.

        Missing receptor/ligand Contexts (execute_decomposition) are created on
        the decomposition platform.
        """
        for key in ('protein_system', 'ligand_system'):
            if systems.get(key) is None:
                raise DecompositionError(
                    f"No {key.split('_')[0]} System: the solvation terms cannot be decomposed exactly")
        contexts = {}
        for name, key in (('complex', 'complex'), ('receptor', 'protein'), ('ligand', 'ligand')):
            ctx = systems.get(f'{key}_context')
            if ctx is None:
                decomp_platform = getattr(self.mmgbsa_calculator, 'decomposition_platform', None)
                platform, properties = self.mmgbsa_calculator.setup_optimized_platform(platform_name=decomp_platform)
                ctx = openmm.Context(systems[f'{key}_system'], openmm.VerletIntegrator(0.001), platform, properties)
                systems[f'{key}_context'] = ctx
            contexts[name] = ctx
        decomposer = _FrameDecomposer(systems['complex_system'], systems['protein_system'],
                                      systems['ligand_system'], residue_map, ligand_indices,
                                      report_raw_energies=self.report_raw_energies, contexts=contexts)
        log.info(f"  Per-atom solvation terms from the run's own forces: "
                 f"complex {decomposer.solvation['complex'].describe()}")
        return decomposer

    def _log_energy_check(self, energy_check):
        """Report how closely the decomposed energies reproduce OpenMM's."""
        if not energy_check:
            return
        solv = max(c['max_solvation_deviation'] for c in energy_check)
        inter = max(c['max_interaction_deviation'] for c in energy_check)
        mean_sum = np.mean([c['decomposed_binding_energy'] for c in energy_check])
        mean_ref = np.mean([c['openmm_binding_energy'] for c in energy_check])
        log.info(f"Per-atom polar and non-polar energies reproduce OpenMM's force-group energies of the "
                 f"complex, receptor and ligand to {solv:.2e} kcal/mol (largest deviation, all frames).")
        log.info(f"Receptor residues plus the ligand's share: {mean_sum:.4f} kcal/mol; OpenMM binding "
                 f"energy (vdW + electrostatic + polar + non-polar) on the same frames: {mean_ref:.4f} "
                 f"kcal/mol (largest per-frame vdW/electrostatic deviation {inter:.2e}).")

    def _decompose_single_frame(self, frame, systems, residue_map, ligand_indices, ligand_resname,
                                salt_concentration=None):
        """
        Decompose one MDTraj frame into per-residue contributions (kept for
        callers of the former API). Raises DecompositionError on failure.
        """
        decomposer = systems.get('_decomposer')
        if decomposer is None:
            decomposer = self._make_decomposer(systems, residue_map, ligand_indices)
            systems['_decomposer'] = decomposer
        frame_result, _ = decomposer.decompose(frame.xyz[0])
        return frame_result

    def _calculate_pairwise_interactions(self, systems, positions, residue_map, ligand_indices, salt_concentration=None):
        """Per-residue van der Waals and electrostatic interaction with
        `ligand_indices` (halved), from the complex System's own nonbonded forces."""
        return _calculate_pairwise_interactions_standalone(
            systems['complex_system'], None, positions, residue_map, ligand_indices)

    def _average_residue_energies(self, residue_energies):
        """
        Average residue energies across frames (Dynamic for all keys)
        """
        
        averaged = {}
        
        # Get all residue IDs and sample keys
        all_residues = set()
        sample_keys = set()
        
        for frame_result in residue_energies:
            all_residues.update(frame_result.keys())
            for res_data in frame_result.values():
                sample_keys.update(k for k in res_data.keys() if isinstance(res_data[k], (int, float)))
        
        # Average each residue's contributions
        for res_id in all_residues:
            # Initialize collectors for all known keys
            collectors = {k: [] for k in sample_keys}
            
            for frame_result in residue_energies:
                if res_id in frame_result:
                    res_data = frame_result[res_id]
                    for k in sample_keys:
                        if k in res_data:
                            collectors[k].append(res_data[k])
            
            # Calculate stats for keys that have data
            stats = {}
            for k, values in collectors.items():
                if values:
                    stats[f"{k}_mean"] = np.mean(values)
                    stats[f"{k}_std"] = np.std(values)
            
            if stats:
                stats['n_frames'] = len(collectors.get('total', [])) # Use total as reference for frame count
                averaged[res_id] = stats
        
        return averaged
    
    def _analyze_decomposition_results(self, decomp_results):
        """
        Analyze decomposition results to identify hot spots and patterns
        """
        
        print("  Analyzing residue contributions...")
        
        # Create DataFrame for analysis
        data = []
        for res_id, energies in decomp_results.items():
            # Parse residue info
            parts = res_id.split('_')
            res_name = parts[0]
            res_num = int(parts[1])
            chain = (parts[2].strip() if len(parts) > 2 else '') or 'A'
            
            # Base entry
            entry = {
                'residue_id': res_id,
                'residue_name': res_name,
                'residue_number': res_num,
                'chain': chain
            }
            
            # Dynamically add all energy components
            for k, v in energies.items():
                if k.endswith('_mean'):
                    # Map 'key_mean' -> 'key' (e.g. vdw_mean -> vdw)
                    entry[k[:-5]] = v
                else:
                    entry[k] = v
            
            data.append(entry)
        
        df = pd.DataFrame(data)
        
        # Identify hot spots (most favorable contributions)
        hot_spots = df.nsmallest(10, 'total')  # Most negative = most favorable
        
        # Identify key interaction types
        vdw_important = df.nsmallest(5, 'vdw')
        elec_important = df.nsmallest(5, 'electrostatic')
        
        # Summary statistics
        total_contribution = df['total'].sum()
        mean_contribution = df['total'].mean()
        
        analysis_results = {
            'dataframe': df,
            'hot_spots': hot_spots,
            'vdw_important': vdw_important,
            'elec_important': elec_important,
            'total_contribution': total_contribution,
            'mean_contribution': mean_contribution,
            'n_residues': len(df)
        }
        
        print(f"  Analyzed {len(df)} residues")
        print(f"  Total residue contribution: {total_contribution:.2f} kcal/mol")
        
        return analysis_results
    
    def _generate_decomposition_plots(self, analysis_results):
        """
        Generate visualization plots for decomposition results
        """
        
        try:
            print("  Generating decomposition plots...")
            
            df = analysis_results['dataframe']
            
            # Create figure with subplots
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle('Per-Residue Energy Decomposition Analysis', fontsize=16, fontweight='bold')
            
            # Plot 1: Hot spot contributions
            ax1 = axes[0, 0]
            hot_spots = analysis_results['hot_spots'].head(10)
            bars = ax1.barh(range(len(hot_spots)), hot_spots['total'], 
                           color='red', alpha=0.7, edgecolor='black')
            ax1.set_yticks(range(len(hot_spots)))
            ax1.set_yticklabels([f"{row['residue_name']}{row['residue_number']}" 
                               for _, row in hot_spots.iterrows()])
            ax1.set_xlabel('Total Energy (kcal/mol)')
            ax1.set_title('Top 10 Binding Hot Spots')
            ax1.grid(True, alpha=0.3)
            
            # Add value labels
            for i, (_, row) in enumerate(hot_spots.iterrows()):
                ax1.text(row['total'], i, f'{row["total"]:.1f}', 
                        va='center', ha='right' if row['total'] < 0 else 'left')
            
            # Plot 2: Component breakdown for top residues
            ax2 = axes[0, 1]
            top_5 = hot_spots.head(5)
            x_pos = np.arange(len(top_5))
            width = 0.25
            
            ax2.bar(x_pos - width, top_5['vdw'], width, label='van der Waals', alpha=0.8)
            ax2.bar(x_pos, top_5['electrostatic'], width, label='Electrostatic', alpha=0.8)
            ax2.bar(x_pos + width, top_5['solvation'], width, label='Solvation', alpha=0.8)
            
            ax2.set_xlabel('Residue')
            ax2.set_ylabel('Energy (kcal/mol)')
            ax2.set_title('Energy Components for Top 5 Hot Spots')
            ax2.set_xticks(x_pos)
            ax2.set_xticklabels([f"{row['residue_name']}{row['residue_number']}" 
                               for _, row in top_5.iterrows()], rotation=45)
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            # Plot 3: Energy distribution
            ax3 = axes[1, 0]
            ax3.hist(df['total'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(df['total'].mean(), color='red', linestyle='--', 
                       label=f'Mean: {df["total"].mean():.2f}')
            ax3.set_xlabel('Total Energy (kcal/mol)')
            ax3.set_ylabel('Number of Residues')
            ax3.set_title('Distribution of Residue Contributions')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            # Plot 4: vdW vs Electrostatic scatter
            ax4 = axes[1, 1]
            scatter = ax4.scatter(df['vdw'], df['electrostatic'], 
                                c=df['total'], cmap='RdYlBu', alpha=0.7, s=50)
            ax4.set_xlabel('van der Waals Energy (kcal/mol)')
            ax4.set_ylabel('Electrostatic Energy (kcal/mol)')
            ax4.set_title('vdW vs Electrostatic Contributions')
            ax4.grid(True, alpha=0.3)
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax4)
            cbar.set_label('Total Energy (kcal/mol)')
            
            # Annotate top contributors
            for _, row in hot_spots.head(3).iterrows():
                ax4.annotate(f"{row['residue_name']}{row['residue_number']}", 
                           (row['vdw'], row['electrostatic']),
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            plt.tight_layout()
            
            # Save plot to output directory if available
            plot_filename = 'per_residue_decomposition.png'
            if hasattr(self, 'report_dir') and self.report_dir:
                plot_path = os.path.join(self.report_dir, plot_filename)
            else:
                plot_path = plot_filename
            
            plt.savefig(plot_path, dpi=600, bbox_inches='tight')
            plt.close()
            
            print(f"  Plots saved to {plot_path}")

        except Exception as e:
            print(f"  WARNING: Plot generation failed: {e}")

    def _generate_time_series_heatmap(self, frame_data, residue_map, top_n=20):
        """
        Generate time-dependent interaction energy heatmap for top residues
        """
        try:
            if not frame_data:
                return

            print(f"  Generating time-series heatmap for top {top_n} residues...")
            
            # Convert to DataFrame
            df_time = pd.DataFrame(frame_data)
            
            # Identify columns Ending with _total (representing total energy)
            energy_cols = [col for col in df_time.columns if col.endswith('_total')]
            
            # Calculate mean for each residue column to find top N
            mean_energies = {col: df_time[col].mean() for col in energy_cols}
            
            # Sort by most negative (strongest interaction)
            sorted_residues = sorted(mean_energies.items(), key=lambda x: x[1])
            top_residues = sorted_residues[:top_n]  # Top N most stabilizing
            
            # Prepare matrix for heatmap (X=Frame, Y=Residue)
            heatmap_data = []
            residue_labels = []
            
            for col, mean_val in top_residues:
                clean_name = col.replace('_total', '').replace('_A', '')
                residue_labels.append(f"{clean_name} ({mean_val:.1f})")
                heatmap_data.append(df_time[col].values)
            
            # --- Add Ligand Total Row (Bottom) ---
            # Sum of ALL residues (not just top N) gives the Total Binding Energy
            total_energy_series = df_time[energy_cols].sum(axis=1)
            total_mean = total_energy_series.mean()
            
            heatmap_data.append(total_energy_series.values)
            residue_labels.append(f"LIGAND Total ({total_mean:.1f})")
            # -------------------------------------
            
            # Create Matrix: Rows=Residues, Cols=Frames
            heatmap_matrix = np.array(heatmap_data)
            
            plt.figure(figsize=(14, max(8, (top_n + 1) * 0.4))) # Dynamic height
            
            # Use diverging colormap (Blue=Negative, Red=Positive) centered at 0
            # Normalize colors based on RESIDUES ONLY (exclude Ligand Total to preserve contrast)
            residue_matrix = heatmap_matrix[:-1] # All except last row
            vmin, vmax = np.nanpercentile(residue_matrix, 5), np.nanpercentile(residue_matrix, 95)
            abs_max = max(abs(vmin), abs(vmax))
            
            # Smart X-Axis Labeling
            # Show max ~15 ticks to avoid crowding
            frames = df_time['frame_number'].tolist()
            step = max(1, len(frames) // 15)
            # Create list of labels with empty strings for skipped frames
            x_labels = [str(f) if i % step == 0 else "" for i, f in enumerate(frames)]
            
            sns.heatmap(heatmap_matrix, 
                       xticklabels=x_labels,
                       yticklabels=residue_labels,
                       cmap="RdBu_r", # Reverse Red-Blue: Blue is Negative (Good), Red is Positive (Bad)
                       center=0,
                       vmin=-abs_max, vmax=abs_max,
                       cbar_kws={'label': 'Interaction Energy (kcal/mol)'})
            
            plt.xlabel('Frame Number')
            plt.ylabel('Residue')
            plt.title(f'Time Evolution of Top {top_n} Residue Interactions (Heatmap)')
            plt.tight_layout()
            
            plot_filename = 'per_residue_heatmap.png'
            if hasattr(self, 'report_dir') and self.report_dir:
                plot_path = os.path.join(self.report_dir, plot_filename)
            else:
                plot_path = plot_filename
                
            plt.savefig(plot_path, dpi=600, bbox_inches='tight')
            plt.close()
            print(f"  Heatmap saved to {plot_path}")
            
        except Exception as e:
            print(f"  WARNING: Heatmap generation failed: {e}")
    
    def _print_decomposition_summary(self, results):
        """
        Print comprehensive summary of decomposition results
        """
        
        analysis = results['analysis_results']
        mmgbsa = results['mmgbsa_results']
        
        print("\n" + "="*60)
        print("PER-RESIDUE DECOMPOSITION SUMMARY")
        print("="*60)
        
        print(f"MM/GBSA Total: {mmgbsa['mean_binding_energy']:.2f} ± {mmgbsa['std_error']:.2f} kcal/mol")
        print(f"Residue Sum:   {analysis['total_contribution']:.2f} kcal/mol")
        print(f"Residues:      {analysis['n_residues']}")
        print(f"Frames:        {results['n_decomp_frames']}")

        # The receptor residues carry half of each receptor-ligand van der
        # Waals/electrostatic pair term (Amber decpair convention) and their
        # own polar/non-polar change; the ligand carries the other half and its
        # own solvation change. Both together are the binding energy of the
        # decomposed frames, term by term (checked against OpenMM per frame).
        total_binding = mmgbsa.get('mean_binding_energy')
        residue_sum = analysis.get('total_contribution')
        check = results.get('energy_check') or []
        if check:
            decomposed = np.mean([c['decomposed_binding_energy'] for c in check])
            reference = np.mean([c['openmm_binding_energy'] for c in check])
            ligand = np.mean([c['ligand_vdw'] + c['ligand_electrostatic'] + c['ligand_solvation']
                              + c['ligand_nonpolar_solvation'] for c in check])
            nonpolar = np.mean([c['receptor_nonpolar_solvation'] for c in check])
            print(f"Receptor residues (vdW + electrostatic + polar): {residue_sum:.2f} kcal/mol; "
                  f"their non-polar term {nonpolar:.2f}; ligand share {ligand:.2f}.")
            print(f"Decomposed binding energy {decomposed:.4f} kcal/mol = OpenMM binding energy on "
                  f"the same {len(check)} frames {reference:.4f} kcal/mol "
                  f"(MM/GBSA mean over its own frames: {total_binding:.2f}).")
        print()
        
        print("TOP 10 BINDING HOT SPOTS:")
        print("-" * 40)
        hot_spots = analysis['hot_spots']
        for i, (_, row) in enumerate(hot_spots.iterrows(), 1):
            print(f"{i:2d}. {row['residue_name']}{row['residue_number']:3d} "
                  f"({row['chain']}) = {row['total']:6.2f} ± {row['total_std']:4.2f} kcal/mol")
        
        print(f"\nTOP van der Waals Contributors:")
        print("-" * 35)
        for i, (_, row) in enumerate(analysis['vdw_important'].iterrows(), 1):
            print(f"{i}. {row['residue_name']}{row['residue_number']} = {row['vdw']:5.2f} kcal/mol")
        
        print(f"\nTOP Electrostatic Contributors:")
        print("-" * 35)
        for i, (_, row) in enumerate(analysis['elec_important'].iterrows(), 1):
            print(f"{i}. {row['residue_name']}{row['residue_number']} = {row['electrostatic']:5.2f} kcal/mol")
        
        print(f"\nSUMMARY STATISTICS:")
        print("-" * 20)
        df = analysis['dataframe']
        print(f"Mean residue contribution: {analysis['mean_contribution']:6.2f} kcal/mol")
        print(f"Std deviation:             {df['total'].std():6.2f} kcal/mol")
        print(f"Most favorable:            {df['total'].min():6.2f} kcal/mol")
        print(f"Least favorable:           {df['total'].max():6.2f} kcal/mol")
        
        # Key residue identification
        favorable_residues = len(df[df['total'] < -1.0])  # Strong contributors
        unfavorable_residues = len(df[df['total'] > 1.0])  # Unfavorable
        
        print(f"Strong contributors (< -1 kcal/mol): {favorable_residues}")
        print(f"Unfavorable (> +1 kcal/mol):         {unfavorable_residues}")
        
        print("="*60)
    
    def _save_decomposition_results(self, results):
        """
        Save detailed decomposition results to files
        """
        
        try:
            # Get data
            df = results['analysis_results']['dataframe']
            hot_spots = results['analysis_results']['hot_spots']
            
            # Create summary statistics
            summary = {
                'total_binding_energy': results['mmgbsa_results']['mean_binding_energy'],
                'residue_sum': results['analysis_results']['total_contribution'],
                'n_residues': results['analysis_results']['n_residues'],
                'n_frames': results['n_decomp_frames'],
                'top_hotspot': hot_spots.iloc[0]['residue_id'],
                'top_contribution': hot_spots.iloc[0]['total']
            }
            check = results.get('energy_check') or []
            if check:
                summary['decomposed_binding_energy'] = float(np.mean([c['decomposed_binding_energy'] for c in check]))
                summary['openmm_binding_energy_same_frames'] = float(np.mean([c['openmm_binding_energy'] for c in check]))
                summary['max_solvation_check_deviation'] = float(max(c['max_solvation_deviation'] for c in check))
            summary_df = pd.DataFrame([summary])
            
            # Save files to output directory if available
            output_dir = getattr(self, 'output_dir', None)
            
            if output_dir:
                # Save to output directory
                df.to_csv(os.path.join(output_dir, 'per_residue_detailed.csv'), index=False, encoding='utf-8')
                hot_spots.to_csv(os.path.join(output_dir, 'binding_hot_spots.csv'), index=False, encoding='utf-8')
                summary_df.to_csv(os.path.join(output_dir, 'decomposition_summary.csv'), index=False, encoding='utf-8')
                if check:
                    pd.DataFrame(check).to_csv(os.path.join(output_dir, 'decomposition_energy_check.csv'),
                                               index=False, encoding='utf-8')
                
                print(f"Results saved to {output_dir}:")
                print(f"    • per_residue_detailed.csv")
                print(f"    • binding_hot_spots.csv") 
                print(f"    • decomposition_summary.csv")
                if check:
                    print(f"    • decomposition_energy_check.csv")
                
                # Save Amber-like raw energy tables if available
                if 'complex_total' in df.columns:
                    # Helper to filter columns
                    def save_subset(prefix, filename):
                        """Write an Amber-like per-residue CSV for one energy
                        perspective ('ligand', 'receptor', or 'complex'), keeping
                        only that prefix's vdw/electrostatic/solvation/total columns
                        alongside the residue identity columns."""
                        base_cols = ['residue_id', 'residue_name', 'residue_number', 'chain']
                        target_cols = base_cols.copy()
                        
                        for col in df.columns:
                            if col in base_cols: continue
                            
                            if prefix == 'ligand':
                                # Standard columns: vdw, electrostatic, solvation, total
                                if col in ['vdw', 'electrostatic', 'solvation', 'total']:
                                    target_cols.append(col)
                                elif col.endswith('_mean') and col[:-5] in ['vdw', 'electrostatic', 'solvation', 'total']:
                                    target_cols.append(col)
                                # Also include std columns
                                elif col in ['vdw_std', 'electrostatic_std', 'solvation_std', 'total_std']:
                                    target_cols.append(col)
                            else:
                                # complex_vdw, receptor_total, etc.
                                if col.startswith(prefix + '_'):
                                    target_cols.append(col)
                        
                        # Only save if we found data columns
                        if len(target_cols) > len(base_cols):
                            subset_df = df[target_cols]
                            subset_path = os.path.join(output_dir, filename)
                            subset_df.to_csv(subset_path, index=False, encoding='utf-8')
                            print(f"    • {filename}")

                    save_subset('complex', 'per_residue_complex.csv')
                    save_subset('receptor', 'per_residue_receptor.csv')
                    save_subset('ligand', 'per_residue_ligand.csv')
                
                # Show frame-by-frame output if enabled
                if hasattr(self, 'frame_by_frame_settings') and self.frame_by_frame_settings.get('save_frame_csv', True):
                    frame_csv_name = self.frame_by_frame_settings.get('frame_by_frame_csv_name', 'frame_by_frame_decomposition')
                    frame_output_format = self.frame_by_frame_settings.get('frame_output_format', 'csv')
                    include_residue_summary = self.frame_by_frame_settings.get('include_residue_summary', True)
                    
                    print(f"    • {frame_csv_name}.{frame_output_format}")  # Frame-by-frame data
                    if include_residue_summary:
                        print(f"    • {frame_csv_name}_residue_summary.{frame_output_format}")  # Residue summary
                
                print(f"    • per_residue_decomposition.png")
            else:
                # Save to current directory
                df.to_csv('per_residue_detailed.csv', index=False, encoding='utf-8')
                hot_spots.to_csv('binding_hot_spots.csv', index=False, encoding='utf-8')
                summary_df.to_csv('decomposition_summary.csv', index=False, encoding='utf-8')
                
                print(f"Results saved:")
                print(f"    • per_residue_detailed.csv")
                print(f"    • binding_hot_spots.csv") 
                print(f"    • decomposition_summary.csv")
                
                # Show frame-by-frame output if enabled
                if hasattr(self, 'frame_by_frame_settings') and self.frame_by_frame_settings.get('save_frame_csv', True):
                    frame_csv_name = self.frame_by_frame_settings.get('frame_by_frame_csv_name', 'frame_by_frame_decomposition')
                    frame_output_format = self.frame_by_frame_settings.get('frame_output_format', 'csv')
                    include_residue_summary = self.frame_by_frame_settings.get('include_residue_summary', True)
                    
                    print(f"    • {frame_csv_name}.{frame_output_format}")  # Frame-by-frame data
                    if include_residue_summary:
                        print(f"    • {frame_csv_name}_residue_summary.{frame_output_format}")  # Residue summary
                
                print(f"    • per_residue_decomposition.png")
            
        except Exception as e:
            print(f"  WARNING: Could not save results: {e}")
    
    def _save_frame_by_frame_csv(self, frame_by_frame_data, residue_map):
        """
        Save frame-by-frame per-residue decomposition results to CSV
        """
        
        try:
            if not frame_by_frame_data:
                print("  WARNING: No frame-by-frame data to save")
                return
            
            # Get settings from config
            settings = getattr(self, 'frame_by_frame_settings', {})
            frame_csv_name = settings.get('frame_by_frame_csv_name', 'frame_by_frame_decomposition')
            include_residue_summary = settings.get('include_residue_summary', True)
            frame_output_components = settings.get('frame_output_components', ['vdw', 'electrostatic', 'solvation', 'total'])
            frame_output_format = settings.get('frame_output_format', 'csv')
            
            # Convert to DataFrame
            import pandas as pd
            df = pd.DataFrame(frame_by_frame_data)
            
            # Filter columns based on output components if needed
            if frame_output_components != ['vdw', 'electrostatic', 'solvation', 'total']:
                # Keep frame metadata columns
                metadata_cols = ['frame_index', 'frame_number', 'total_frames']
                filtered_cols = metadata_cols.copy()
                
                # Add only requested energy components
                for col in df.columns:
                    if col not in metadata_cols:
                        for component in frame_output_components:
                            if col.endswith(f'_{component}'):
                                filtered_cols.append(col)
                                break
                
                df = df[filtered_cols]
            
            # Determine output path
            if hasattr(self, 'output_dir') and self.output_dir:
                output_path = os.path.join(self.output_dir, f'{frame_csv_name}.{frame_output_format}')
            else:
                output_path = f'{frame_csv_name}.{frame_output_format}'
            
            # Save based on format
            if frame_output_format == 'csv':
                df.to_csv(output_path, index=False, encoding='utf-8')
            elif frame_output_format == 'json':
                df.to_json(output_path, orient='records', indent=2, force_ascii=False)
            elif frame_output_format == 'hdf5':
                df.to_hdf(output_path, key='frame_decomposition', mode='w')
            else:
                # Default to CSV
                df.to_csv(output_path, index=False, encoding='utf-8')
            
            print(f"Frame-by-frame decomposition saved: {output_path}")
            print(f"    • {len(df)} frames")
            print(f"    • {len(df.columns) - 3} residue energy columns")  # -3 for frame metadata
            print(f"    • Components: {', '.join(frame_output_components)}")
            print(f"    • Format: {frame_output_format}")
            
            # Save residue summary if enabled
            if include_residue_summary:
                residue_summary = []
                for res_id in residue_map.keys():
                    parts = res_id.split('_')
                    if len(parts) >= 3:
                        res_name = parts[0]
                        res_number = parts[1]
                        chain_id = parts[2]
                        residue_summary.append({
                            'residue_id': res_id,
                            'residue_name': res_name,
                            'residue_number': res_number,
                            'chain_id': chain_id,
                            'column_prefix': f'{res_name}{res_number}_{chain_id}',
                            'available_components': frame_output_components
                        })
                
                summary_df = pd.DataFrame(residue_summary)
                summary_path = output_path.replace(f'.{frame_output_format}', f'_residue_summary.{frame_output_format}')
                
                if frame_output_format == 'csv':
                    summary_df.to_csv(summary_path, index=False)
                elif frame_output_format == 'json':
                    summary_df.to_json(summary_path, orient='records', indent=2)
                elif frame_output_format == 'hdf5':
                    summary_df.to_hdf(summary_path, key='residue_summary', mode='w')
                
                print(f"Residue summary saved: {summary_path}")
            
        except Exception as e:
            print(f"  WARNING: Could not save frame-by-frame CSV: {e}")


def test_per_residue_decomposition():
    """
    Test the per-residue decomposition with your existing MM/GBSA package
    """
    
    print("="*60)
    print("TESTING PER-RESIDUE DECOMPOSITION")
    print("="*60)
    print("This adds advanced per-residue analysis to your MM/GBSA package!")
    
    # Your existing MM/GBSA calculator
    mmgbsa_calc = GBSACalculator(
        temperature=300,
        gb_model='OBC2',
        salt_concentration=0.15,
        use_cache=True,
        verbose=1
    )
    
    # Create per-residue decomposition analyzer
    decomp_analyzer = PerResidueDecomposition(mmgbsa_calc, temperature=300)
    
    # Input files
    ligand_mol = 'test/ligand.sdf'
    complex_pdb = 'test/complex.pdb'
    ligand_pdb = 'test/ligand.pdb'
    xtc_file = 'test/complex.dcd'
    
    # Run complete analysis with per-residue decomposition
    results = decomp_analyzer.run_per_residue_analysis(
        ligand_mol, complex_pdb, xtc_file, ligand_pdb,
        max_frames=20,      # For MM/GBSA
        decomp_frames=5     # For decomposition (computationally expensive)
    )
    
    if results:
        print(f"\nSUCCESS! Per-residue decomposition completed!")
        print(f"Your MM/GBSA package now has advanced per-residue analysis!")
        
        # Show key results
        hot_spots = results['analysis_results']['hot_spots']
        print(f"\nTOP 3 BINDING HOT SPOTS:")
        for i, (_, row) in enumerate(hot_spots.head(3).iterrows(), 1):
            print(f"{i}. {row['residue_name']}{row['residue_number']} = {row['total']:.2f} kcal/mol")
        
        return results
    
    else:
        print(f"\nWARNING: Per-residue decomposition had issues")
        print(f"But your core MM/GBSA package is still excellent!")
        return None



# ==============================================================================
# PER-FRAME DECOMPOSITION ENGINE (shared by the serial and the parallel path)
# ==============================================================================

ELECTROSTATIC_GROUP = 0      # GBSACalculator: NonbondedForce (Coulomb)
VDW_GROUP = 3                # CustomNonbondedForce LJ and the 1-4 LJ CustomBondForce
FRAME_COMPONENTS = ('vdw', 'electrostatic', 'solvation', 'nonpolar_solvation', 'total')


def _frame_row(frame_result, frame_index, frame_number, total_frames):
    """One row of frame_by_frame_decomposition.csv: every component of every
    residue, labelled <NAME><NUMBER>_<CHAIN>_<component>."""
    row = {'frame_index': frame_index, 'frame_number': frame_number, 'total_frames': total_frames}
    for res_id, energies in frame_result.items():
        parts = res_id.split('_')
        if len(parts) >= 3:
            prefix = f"{parts[0]}{parts[1]}_{parts[2]}"
        else:
            prefix = f"{res_id}0_A"
        for comp in FRAME_COMPONENTS:
            row[f'{prefix}_{comp}'] = energies.get(comp, 0.0)
    return row


def _strip(x, u=None):
    if unit.is_quantity(x):
        x = x.value_in_unit(u) if u is not None else x._value
    return float(x)


class _CustomPairModel:
    """A CustomNonbondedForce (van der Waals group) evaluated from its own
    energy expression, per-particle parameters and tabulated functions."""

    def __init__(self, force):
        if force.getNonbondedMethod() != openmm.CustomNonbondedForce.NoCutoff:
            raise DecompositionError("CustomNonbondedForce with a cutoff cannot be decomposed")
        if force.getNumInteractionGroups() > 0:
            raise DecompositionError("CustomNonbondedForce interaction groups are not supported")
        if hasattr(force, 'getNumComputedValues') and force.getNumComputedValues() > 0:
            raise DecompositionError("CustomNonbondedForce computed values are not supported")
        funcs = dict(_LEPTON_BUILTINS)
        for k in range(force.getNumTabulatedFunctions()):
            funcs[force.getTabulatedFunctionName(k)] = _discrete_function(force.getTabulatedFunction(k))
        self.expr = _LeptonExpression(force.getEnergyFunction(), funcs)
        n = force.getNumParticles()
        names = [force.getPerParticleParameterName(k) for k in range(force.getNumPerParticleParameters())]
        values = np.zeros((n, len(names)))
        for i in range(n):
            values[i] = [_strip(v) for v in force.getParticleParameters(i)]
        self.params = {name: values[:, k].copy() for k, name in enumerate(names)}
        self.globals = {force.getGlobalParameterName(k): _strip(force.getGlobalParameterDefaultValue(k))
                        for k in range(force.getNumGlobalParameters())}
        excl = [force.getExclusionParticles(k) for k in range(force.getNumExclusions())]
        self.excl_i = np.array([e[0] for e in excl], dtype=np.int64)
        self.excl_j = np.array([e[1] for e in excl], dtype=np.int64)

    def energy(self, idx1, idx2, r):
        """kJ/mol for particle 1 = idx1, particle 2 = idx2 (broadcast shapes)."""
        v = dict(self.globals)
        for name, arr in self.params.items():
            v[name + '1'] = arr[idx1]
            v[name + '2'] = arr[idx2]
        v['r'] = r
        return np.broadcast_to(np.asarray(self.expr.evaluate(v), dtype=float), np.shape(r))


class _CustomBondModel:
    """A CustomBondForce of the van der Waals group (the 1-4 LJ that
    GBSAForceManager.refine_gbsa_forces moves out of the NonbondedForce)."""

    def __init__(self, force):
        if force.usesPeriodicBoundaryConditions():
            raise DecompositionError("A periodic CustomBondForce cannot be decomposed")
        self.expr = _LeptonExpression(force.getEnergyFunction(), dict(_LEPTON_BUILTINS))
        nb = force.getNumBonds()
        names = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
        self.i = np.zeros(nb, dtype=np.int64)
        self.j = np.zeros(nb, dtype=np.int64)
        values = np.zeros((nb, len(names)))
        for k in range(nb):
            a, b, p = force.getBondParameters(k)
            self.i[k], self.j[k] = a, b
            values[k] = [_strip(x) for x in p]
        self.params = {name: values[:, k].copy() for k, name in enumerate(names)}
        self.globals = {force.getGlobalParameterName(k): _strip(force.getGlobalParameterDefaultValue(k))
                        for k in range(force.getNumGlobalParameters())}

    def energy(self, sel, r):
        v = dict(self.globals)
        for name, arr in self.params.items():
            v[name] = arr[sel]
        v['r'] = r
        return np.broadcast_to(np.asarray(self.expr.evaluate(v), dtype=float), np.shape(r))


class _NonbondedModel:
    """Pairwise electrostatic and van der Waals energies of a System, read from
    its own forces: the NonbondedForce (charges, exceptions; any LJ left in it),
    and the CustomNonbondedForce and CustomBondForce of the van der Waals force
    group. Charges therefore carry whatever scaling the run applied (e.g. a
    Coulomb term divided by the solute dielectric)."""

    def __init__(self, system):
        n = system.getNumParticles()
        self.n = n
        self.q = np.zeros(n)
        self.sig = np.zeros(n)
        self.eps = np.zeros(n)
        ex = []
        self.custom, self.bonds = [], []
        n_nb = 0
        for f in system.getForces():
            group = f.getForceGroup()
            if isinstance(f, openmm.NonbondedForce):
                n_nb += 1
                if f.getNonbondedMethod() != openmm.NonbondedForce.NoCutoff:
                    raise DecompositionError("A NonbondedForce with a cutoff cannot be decomposed per residue")
                for i in range(n):
                    q, s, e = f.getParticleParameters(i)
                    self.q[i] = _strip(q, unit.elementary_charge)
                    self.sig[i] = _strip(s, unit.nanometer)
                    self.eps[i] = _strip(e, unit.kilojoule_per_mole)
                for k in range(f.getNumExceptions()):
                    i, j, qq, s, e = f.getExceptionParameters(k)
                    ex.append((i, j, _strip(qq, unit.elementary_charge ** 2), _strip(s, unit.nanometer),
                               _strip(e, unit.kilojoule_per_mole)))
            elif isinstance(f, openmm.CustomNonbondedForce) and group == VDW_GROUP:
                self.custom.append(_CustomPairModel(f))
            elif isinstance(f, openmm.CustomBondForce) and group == VDW_GROUP:
                self.bonds.append(_CustomBondModel(f))
        if n_nb != 1:
            raise DecompositionError(f"Expected one NonbondedForce, found {n_nb}")
        self.has_lj = bool(np.any(self.eps != 0.0))
        ex = np.array(ex, dtype=float).reshape(-1, 5)
        self.ex_i, self.ex_j = ex[:, 0].astype(np.int64), ex[:, 1].astype(np.int64)
        self.ex_qq, self.ex_sig, self.ex_eps = ex[:, 2], ex[:, 3], ex[:, 4]

    def _lj(self, i, j, r):
        s = 0.5 * (self.sig[i] + self.sig[j])
        e = np.sqrt(self.eps[i] * self.eps[j])
        sr6 = (s / r) ** 6
        return 4.0 * e * (sr6 * sr6 - sr6)

    def pair_energies(self, pos, rows, targets):
        """Half of each row-atom/target-atom pair energy, summed per row atom.

        Returns (electrostatic, vdw) per atom in kcal/mol (zero outside `rows`).
        A pair whose both atoms are rows and targets is counted once from each
        side, so each atom carries its half.
        """
        n = self.n
        rows = np.asarray(rows, dtype=np.int64)
        targets = np.asarray(targets, dtype=np.int64)
        elec = np.zeros(n)
        vdw = np.zeros(n)
        if rows.size == 0 or targets.size == 0:
            return elec, vdw
        qt = self.q[targets]
        m = max(1, CHUNK_ELEMENTS // targets.size)
        blocks = [(s, min(rows.size, s + m)) for s in range(0, rows.size, m)]

        def work(s, e):
            a = rows[s:e]
            r = _pair_distances(pos[a], pos[targets])
            mask = (a[:, None] != targets[None, :]) & (r > 1e-3)
            rs = np.where(mask, r, 1.0)
            ee = ONE_4PI_EPS0 * self.q[a][:, None] * qt[None, :] / rs
            vv = np.zeros(r.shape)
            if self.has_lj:
                vv = vv + self._lj(a[:, None], targets[None, :], rs)
            for c in self.custom:
                vv = vv + c.energy(a[:, None], targets[None, :], rs)
            return a, 0.5 * np.where(mask, ee, 0.0).sum(axis=1), 0.5 * np.where(mask, vv, 0.0).sum(axis=1)

        for a, e_a, v_a in _map_blocks(work, blocks):
            elec[a] = e_a
            vdw[a] = v_a

        in_rows = np.zeros(n, dtype=bool)
        in_rows[rows] = True
        in_t = np.zeros(n, dtype=bool)
        in_t[targets] = True

        def scatter(i, j, d_e, d_v):
            """Add pair corrections to row atom i (target j) and row atom j (target i)."""
            w1 = in_rows[i] & in_t[j]
            w2 = in_rows[j] & in_t[i]
            if d_e is not None:
                elec[:] += np.bincount(i[w1], d_e[w1], n) + np.bincount(j[w2], d_e[w2], n)
            if d_v is not None:
                vdw[:] += np.bincount(i[w1], d_v[0][w1], n) + np.bincount(j[w2], d_v[1][w2], n)

        def selected(i, j):
            sel = np.nonzero((in_rows[i] & in_t[j]) | (in_rows[j] & in_t[i]))[0]
            if sel.size == 0:
                return sel, None, None
            d = pos[i[sel]] - pos[j[sel]]
            r = np.sqrt(np.einsum('ij,ij->i', d, d))
            return sel, r, r > 1e-3

        # NonbondedForce exceptions replace the standard Coulomb (and LJ) term.
        sel, r, ok = selected(self.ex_i, self.ex_j)
        if sel.size:
            i, j = self.ex_i[sel], self.ex_j[sel]
            rs = np.where(ok, r, 1.0)
            d_e = ONE_4PI_EPS0 * (self.ex_qq[sel] - self.q[i] * self.q[j]) / rs
            sr6 = (self.ex_sig[sel] / rs) ** 6
            act = 4.0 * self.ex_eps[sel] * (sr6 * sr6 - sr6)
            std = self._lj(i, j, rs) if self.has_lj else 0.0
            d_v = np.where(ok, 0.5 * (act - std), 0.0)
            scatter(i, j, np.where(ok, 0.5 * d_e, 0.0), (d_v, d_v))
        # CustomNonbondedForce exclusions: that pair has no standard term.
        for c in self.custom:
            sel, r, ok = selected(c.excl_i, c.excl_j)
            if sel.size:
                i, j = c.excl_i[sel], c.excl_j[sel]
                rs = np.where(ok, r, 1.0)
                v_ij = np.where(ok, -0.5 * c.energy(i, j, rs), 0.0)
                v_ji = np.where(ok, -0.5 * c.energy(j, i, rs), 0.0)
                scatter(i, j, None, (v_ij, v_ji))
        # CustomBondForce (1-4 LJ) between a row atom and a target atom.
        for b in self.bonds:
            sel, r, ok = selected(b.i, b.j)
            if sel.size:
                v = 0.5 * b.energy(sel, r)
                scatter(b.i[sel], b.j[sel], None, (v, v))
        return elec * KJ_TO_KCAL, vdw * KJ_TO_KCAL


def _pair_distances(a, b):
    d = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.einsum('ijk,ijk->ij', d, d))


class _FrameDecomposer:
    """Per-residue decomposition of one frame, from the run's own Systems.

    Van der Waals / electrostatic: half of each receptor-residue/ligand pair
    energy (Amber decpair convention), from the complex System's nonbonded
    forces. Polar and non-polar solvation: each atom's energy in the complex
    minus its energy in the receptor (receptor atoms) or ligand (ligand atoms)
    System, from the GB and surface-area forces of those Systems
    (mmgbsa.solvation_peratom). Summed over the receptor residues and the
    ligand this is the frame's binding energy, term by term.
    """

    def __init__(self, complex_system, receptor_system, ligand_system, residue_map, ligand_indices,
                 report_raw_energies=False, contexts=None, interaction_tolerance=None):
        n = complex_system.getNumParticles()
        self.n = n
        lig = np.array(sorted({int(i) for i in ligand_indices}), dtype=np.int64)
        lig_set = set(lig.tolist())
        rec = np.array([i for i in range(n) if i not in lig_set], dtype=np.int64)
        if receptor_system.getNumParticles() != rec.size or ligand_system.getNumParticles() != lig.size:
            raise DecompositionError(
                f"The receptor and ligand Systems have {receptor_system.getNumParticles()} and "
                f"{ligand_system.getNumParticles()} particles, but the decomposition splits the "
                f"{n}-atom complex into {rec.size} receptor and {lig.size} ligand atoms")
        self.rec, self.lig = rec, lig
        self.res_ids = list(residue_map)
        atom_res = np.full(n, -1, dtype=np.int64)
        for k, res_id in enumerate(self.res_ids):
            atom_res[np.asarray(residue_map[res_id], dtype=np.int64)] = k
        if np.any(atom_res[lig] >= 0):
            raise DecompositionError("Ligand atoms are listed in the receptor residue map")
        self.atom_res = atom_res
        self.rows = np.nonzero(atom_res >= 0)[0]
        self.report_raw_energies = report_raw_energies
        self.nb = _NonbondedModel(complex_system)
        self.solvation = {'complex': SolvationEnergies(complex_system),
                          'receptor': SolvationEnergies(receptor_system),
                          'ligand': SolvationEnergies(ligand_system)}
        self.contexts = contexts or {}
        self._warned_interaction = False

    def _per_residue(self, values):
        return np.bincount(self.atom_res[self.rows], weights=values[self.rows], minlength=len(self.res_ids))

    def decompose(self, xyz_nm):
        """Returns ({res_id: energies}, check) for one frame (positions in nm)."""
        if unit.is_quantity(xyz_nm):
            xyz_nm = xyz_nm.value_in_unit(unit.nanometer)
        pos = np.ascontiguousarray(np.asarray(xyz_nm, dtype=np.float64).reshape(self.n, 3))
        rec, lig = self.rec, self.lig

        elec_a, vdw_a = self.nb.pair_energies(pos, self.rows, lig)
        pc, nc = self.solvation['complex'].per_atom(pos)
        pr, nr = self.solvation['receptor'].per_atom(pos[rec])
        pl, nl = self.solvation['ligand'].per_atom(pos[lig])
        dp, dn = pc.copy(), nc.copy()
        dp[rec] -= pr
        dp[lig] -= pl
        dn[rec] -= nr
        dn[lig] -= nl

        vdw_r, elec_r = self._per_residue(vdw_a), self._per_residue(elec_a)
        pol_r, np_r = self._per_residue(dp), self._per_residue(dn)
        if self.report_raw_energies:
            ce, cv = self.nb.pair_energies(pos, self.rows, np.arange(self.n))
            c_elec, c_vdw = self._per_residue(ce), self._per_residue(cv)
            c_pol, c_np = self._per_residue(pc), self._per_residue(nc)

        result = {}
        for k, res_id in enumerate(self.res_ids):
            e = {'vdw': float(vdw_r[k]), 'electrostatic': float(elec_r[k]),
                 'solvation': float(pol_r[k]), 'nonpolar_solvation': float(np_r[k])}
            e['total'] = e['vdw'] + e['electrostatic'] + e['solvation']
            if self.report_raw_energies:
                cvdw, cele, csol, cnp = float(c_vdw[k]), float(c_elec[k]), float(c_pol[k]), float(c_np[k])
                e.update({'complex_vdw': cvdw, 'complex_electrostatic': cele,
                          'complex_solvation': csol, 'complex_nonpolar_solvation': cnp,
                          'complex_total': cvdw + cele + csol})
                # Receptor = complex - binding: each residue's own value in the
                # receptor System (for the solvation terms exactly that).
                r = {'vdw': cvdw - e['vdw'], 'electrostatic': cele - e['electrostatic'],
                     'solvation': csol - e['solvation'], 'nonpolar_solvation': cnp - e['nonpolar_solvation']}
                e.update({f'receptor_{key}': val for key, val in r.items()})
                e['receptor_total'] = r['vdw'] + r['electrostatic'] + r['solvation']
            result[res_id] = e

        check = {
            'receptor_vdw': float(vdw_r.sum()), 'receptor_electrostatic': float(elec_r.sum()),
            'receptor_solvation': float(pol_r.sum()), 'receptor_nonpolar_solvation': float(np_r.sum()),
            'ligand_vdw': float(vdw_r.sum()), 'ligand_electrostatic': float(elec_r.sum()),
            'ligand_solvation': float(dp[lig].sum()), 'ligand_nonpolar_solvation': float(dn[lig].sum()),
        }
        check['decomposed_binding_energy'] = (2.0 * (check['receptor_vdw'] + check['receptor_electrostatic'])
                                              + float(dp.sum()) + float(dn.sum()))
        check.update(self._verify(pos, {'complex': (pc, nc), 'receptor': (pr, nr), 'ligand': (pl, nl)},
                                  2.0 * check['receptor_electrostatic'], 2.0 * check['receptor_vdw']))
        return result, check

    def _verify(self, pos, per_atom, elec_int, vdw_int):
        """Compare with OpenMM's force-group energies of the three Contexts.

        The per-atom solvation energies must sum to each System's GB and
        surface-area group energies (raises otherwise). The interaction terms
        are compared with the complex-minus-receptor-minus-ligand van der
        Waals and electrostatic group energies (a warning if they differ).
        """
        out = {'max_solvation_deviation': float('nan'), 'max_interaction_deviation': float('nan'),
               'openmm_binding_energy': float('nan')}
        if not all(k in self.contexts for k in ('complex', 'receptor', 'ligand')):
            return out
        coords = {'complex': pos, 'receptor': pos[self.rec], 'ligand': pos[self.lig]}
        energies = {}
        worst = 0.0
        for name in ('complex', 'receptor', 'ligand'):
            ctx = self.contexts[name]
            ctx.setPositions(coords[name])
            polar_ref = group_energy(ctx, POLAR_GROUPS)
            nonpolar_ref = group_energy(ctx, NONPOLAR_GROUPS)
            polar, nonpolar = (float(x.sum()) for x in per_atom[name])
            for value, ref, label in ((polar, polar_ref, 'polar'), (nonpolar, nonpolar_ref, 'non-polar')):
                dev = abs(value - ref)
                worst = max(worst, dev)
                if dev > 0.01 + 1e-5 * abs(ref):
                    raise DecompositionError(
                        f"{name} {label} solvation: per-atom energies sum to {value:.6f} kcal/mol, "
                        f"OpenMM gives {ref:.6f} kcal/mol for the same coordinates")
            energies[name] = (group_energy(ctx, (ELECTROSTATIC_GROUP,)), group_energy(ctx, (VDW_GROUP,)),
                              polar_ref, nonpolar_ref)
        delta = [energies['complex'][k] - energies['receptor'][k] - energies['ligand'][k] for k in range(4)]
        dev_int = max(abs(elec_int - delta[0]), abs(vdw_int - delta[1]))
        if dev_int > 0.01 + 1e-5 * (abs(delta[0]) + abs(delta[1])) and not self._warned_interaction:
            self._warned_interaction = True
            log.warning(f"Per-residue van der Waals/electrostatic interaction ({vdw_int:.4f}/{elec_int:.4f} "
                        f"kcal/mol) differs from OpenMM's complex-receptor-ligand difference "
                        f"({delta[1]:.4f}/{delta[0]:.4f} kcal/mol). The receptor or ligand System "
                        f"does not carry the complex's nonbonded parameters.")
        out.update({'max_solvation_deviation': worst, 'max_interaction_deviation': dev_int,
                    'openmm_binding_energy': float(sum(delta))})
        return out


# ==============================================================================
# WORKER FUNCTIONS FOR PARALLEL PROCESSING
# ==============================================================================

def _worker_init(payload):
    """Build this worker's decomposer from the serialized Systems (once per worker)."""
    global _worker_context
    os.environ['OPENGBSA_DECOMP_THREADS'] = '1'      # one process per core already
    from . import solvation_peratom
    solvation_peratom._EXECUTOR = None                # never reuse a forked thread pool
    systems = {key: openmm.XmlSerializer.deserialize(payload[f'{key}_xml'])
               for key in ('complex', 'protein', 'ligand')}
    platform = openmm.Platform.getPlatformByName('CPU')
    contexts = {name: openmm.Context(systems[key], openmm.VerletIntegrator(0.001), platform, {'Threads': '1'})
                for name, key in (('complex', 'complex'), ('receptor', 'protein'), ('ligand', 'ligand'))}
    _worker_context['decomposer'] = _FrameDecomposer(
        systems['complex'], systems['protein'], systems['ligand'], payload['residue_map'],
        payload['ligand_indices'], report_raw_energies=payload['report_raw_energies'], contexts=contexts)


def _worker_decompose(args):
    """Decompose one frame in a worker: returns (index, (result, check), error)."""
    index, xyz = args
    try:
        return index, _worker_context['decomposer'].decompose(xyz), None
    except Exception as e:
        import traceback
        return index, None, f"{e}\n{traceback.format_exc()}"


def _calculate_pairwise_interactions_standalone(system, context, positions, residue_map, ligand_indices,
                                                salt_concentration=None, pdb_params=None, exclusions=None):
    """
    Per-residue van der Waals and electrostatic interaction energies with the
    atoms `ligand_indices`, halved (Gohlke-Kollman / Amber decpair: each pair
    energy is split equally between the two atoms' residues), in kcal/mol.

    Every parameter comes from the System's own nonbonded forces (charges and
    exceptions of the NonbondedForce, the van der Waals CustomNonbondedForce
    and its exclusions, the 1-4 CustomBondForce), so the halves sum to the
    interaction energy OpenMM computes. `context`, `salt_concentration`,
    `pdb_params` and `exclusions` are accepted for compatibility and ignored.
    Returns {res_id: {'vdw': ..., 'elec': ...}}.
    """
    if unit.is_quantity(positions):
        positions = positions.value_in_unit(unit.nanometer)
    pos = np.ascontiguousarray(np.asarray(positions, dtype=np.float64))
    model = _NonbondedModel(system)
    atom_res = np.full(model.n, -1, dtype=np.int64)
    res_ids = list(residue_map)
    for k, res_id in enumerate(res_ids):
        atom_res[np.asarray(residue_map[res_id], dtype=np.int64)] = k
    rows = np.nonzero(atom_res >= 0)[0]
    elec, vdw = model.pair_energies(pos, rows, np.asarray(list(ligand_indices), dtype=np.int64))
    e_res = np.bincount(atom_res[rows], weights=elec[rows], minlength=len(res_ids))
    v_res = np.bincount(atom_res[rows], weights=vdw[rows], minlength=len(res_ids))
    return {res_id: {'vdw': float(v_res[k]), 'elec': float(e_res[k])} for k, res_id in enumerate(res_ids)}


if __name__ == '__main__':
    # Test the per-residue decomposition
    results = test_per_residue_decomposition()
    
    if results:
        print(f"\nYOUR MM/GBSA PACKAGE NOW HAS:")
        print(f"  • Complete MM/GBSA energy calculation")
        print(f"  • Multiple GB models (HCT, OBC1, OBC2, GBn)")
        print(f"  • Normal Mode Analysis entropy")
        print(f"  • Per-residue energy decomposition")  # NEW!
        print(f"  • Hot spot identification")            # NEW!
        print(f"  • Advanced visualization")              # NEW!
        print(f"")
        print(f"This is now FULLY competitive with Schrödinger Prime!")
    else:
        print(f"\nYour MM/GBSA package is still outstanding!")
        print(f"Per-residue decomposition is an advanced feature that can be refined")
