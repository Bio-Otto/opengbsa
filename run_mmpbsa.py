#!/usr/bin/env python3
import sys
import yaml
from pathlib import Path
from mmgbsa.core import GBSACalculator
from mmgbsa.core import GBSACalculator

class OutputLogger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a", buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

# Load config
if len(sys.argv) > 1:
    config_path = sys.argv[1]
else:
    config_path = 'test/5IY2/config_universal.yaml'
with open(config_path) as f:
    config = yaml.safe_load(f)

# Validate and Standardize Inputs
try:
    from mmgbsa.inputs import InputManager
    config = InputManager.validate_inputs(config, base_dir=Path(config_path).parent)
except ImportError:
    pass

print(f"Running Bare Metal Analysis with config: {config_path}")

inp = config['input']
sett = config.get('analysis_settings', {})
print(f"DEBUG INPUT: {inp}")
print(f"DEBUG SETTINGS: {sett}")

calc = GBSACalculator(
    temperature=sett.get('temperature', 300),
    verbose=sett.get('verbose', 1),
    gb_model=sett.get('gb_model', 'OBC2'),
    salt_concentration=sett.get('salt_concentration', 0.15),
    charge_method=sett.get('charge_method', 'am1bcc'),
    solute_dielectric=sett.get('solute_dielectric', 1.0),
    solvent_dielectric=sett.get('solvent_dielectric', 78.5),
    entropy_method=sett.get('entropy_method', 'interaction'),
    decomposition_method=sett.get('decomposition_method', 'full'),
    protein_forcefield=config['params'].get('protein_forcefield', 'amber'),
    use_cache=sett.get('use_cache', True),
    visualization_settings=sett.get('visualization', {})
)
# Set Ligand Forcefield (gaff/openff)
calc.set_ligand_forcefield(config.get('params', {}).get('ligand_forcefield', 'openff'))

inp = config['input']

# Auto-Convert Gromacs if needed
topology_path = str(inp['topology'])
if topology_path.endswith('.top'):
    try:
        from mmgbsa.conversion import GromacsPreprocessor
        print(f"ℹ️  Gromacs Native Input detected: {topology_path}")
        print("ℹ️  Initiating automatic conversion to Amber format...")
        
        top_path = Path(topology_path)
        # Heuristic: verify if trajectory is a coordinate file or if a .gro exists
        coord_file = None
        possible_gro = top_path.parent / (top_path.stem + ".gro")
        
        # If trajectory is .gro, use it
        traj_path = str(inp['trajectory'])
        if traj_path.endswith('.gro') or traj_path.endswith('.pdb'):
            coord_file = traj_path
        elif possible_gro.exists():
            coord_file = str(possible_gro)
        else:
             # Search directory for any .gro
             gros = list(top_path.parent.glob("*.gro"))
             if gros:
                 coord_file = str(gros[0])
        
        if not coord_file:
            print("⚠️  Could not find .gro file for conversion. Conversion might fail.")
            coord_file = traj_path # Last resort
            
        print(f"ℹ️  Using coordinates from: {coord_file}")
        
        converted = GromacsPreprocessor.convert_to_amber(topology_path, coord_file)
        
        # Update config with converted files
        inp['topology'] = converted['topology'] # complex_dry
        inp['solvated_topology'] = converted['solvated_topology']
        inp['receptor_topology'] = converted['receptor_topology']
        inp['ligand_topology'] = converted['ligand_topology']
        
        print("✅ Auto-conversion successful. Proceeding with Native Amber Mode.")
        
    except Exception as e:
        print(f"⚠️  Auto-conversion failed: {e}")
        print("   Proceeding with original input (likely to fail if not handled by generic loader)")
# Universal Mapping
topology = inp['topology']
trajectory = inp['trajectory']
receptor_top = inp.get('receptor_topology')
ligand_top = inp.get('ligand_topology')
solvated_top = inp.get('solvated_topology')
ligand_mol = inp.get('ligand_mol')
ligand_pdb = inp.get('ligand_pdb')
max_frames = sett.get('max_frames', 10)
energy_decomp = sett.get('run_per_residue_decomposition', False)
output_dir = Path(config['params']['output_directory'])
output_dir.mkdir(parents=True, exist_ok=True)

# Redirect Output to Log File
sys.stdout = OutputLogger(output_dir / "analysis.log")
sys.stderr = sys.stdout # Redirect stderr to same log
print(f"📝 Logging to: {output_dir / 'analysis.log'}")

print(f"Calling calc.run directly with complex_pdb={topology}")

try:
    results = calc.run(
        ligand_mol=ligand_mol,
        complex_pdb=topology,
        xtc_file=trajectory,
        ligand_pdb=ligand_pdb,
        receptor_topology=receptor_top,
        ligand_topology=ligand_top,
        solvated_topology=solvated_top,
        max_frames=max_frames,
        energy_decomposition=energy_decomp,
        output_dir=output_dir
    )
    if results:
        print("\nAnalysis SUCCESS!")
        print(f"Delta G: {results.get('delta_g', 'N/A')} kcal/mol")
        print(f"Results saved to: {results.get('output_file', 'Unknown')}")
    else:
        print("\nAnalysis returned None.")
except Exception as e:
    print(f"\nAnalysis FAILED: {e}")
    import traceback
    traceback.print_exc()
