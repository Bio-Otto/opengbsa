"""
Configuration Management Module for OpenGBSA

This module handles configuration file loading, validation, and management.

`OPTIONS` is the single table of every configuration section and option the
analysis reads: its default, the values it accepts and its text in the
template that ``opengbsa --create-complete-config`` writes. Validation, the
defaults filled into a run's ``analysis_config.yaml`` and the template are all
derived from it. ``DEPRECATED_ALIASES`` lists the old spellings that still
work, ``REMOVED_SECTIONS``/``REMOVED_OPTIONS`` the options that no longer do
anything; everything else is an unknown option and gives a warning.
"""

import copy
import difflib
import os
import yaml
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Value checks. Each returns the normalised value or raises ValueError with a
# text that says what is allowed.
# --------------------------------------------------------------------------

def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _or_null(text, null):
    return f"{text}, or null ({null})" if null else text


def boolean(allow_int=False):
    """true or false (with allow_int also an integer >= 0: 0 is false, others true)."""
    def check(value):
        if isinstance(value, bool):
            return value
        if allow_int and _is_int(value) and value >= 0:
            return bool(value)
        raise ValueError("must be true or false" + (" (or 1/0)" if allow_int else ""))
    return check


def integer(minimum=None, maximum=None, null=None, also=()):
    """An integer in [minimum, maximum], or one of `also`; null if `null` says what it means."""
    if minimum is not None and maximum is not None:
        text = f"an integer from {minimum} to {maximum}"
    elif minimum is not None:
        text = f"an integer >= {minimum}"
    else:
        text = "an integer"
    if also:
        text = " or ".join(str(a) for a in also) + f" or {text}"
    text = _or_null(text, null)

    def check(value):
        if value is None and null:
            return None
        if _is_int(value) and (value in also or ((minimum is None or value >= minimum)
                                                 and (maximum is None or value <= maximum))):
            return value
        raise ValueError(f"must be {text}")
    return check


def number(greater_than=None, minimum=None, maximum=None, null=None):
    """A number (> greater_than, >= minimum, <= maximum); returned as a float."""
    parts = []
    if greater_than is not None:
        parts.append(f"> {greater_than:g}")
    if minimum is not None:
        parts.append(f">= {minimum:g}")
    if maximum is not None:
        parts.append(f"<= {maximum:g}")
    text = _or_null("a number" + (" " + " and ".join(parts) if parts else ""), null)

    def check(value):
        if value is None and null:
            return None
        if (_is_number(value)
                and (greater_than is None or value > greater_than)
                and (minimum is None or value >= minimum)
                and (maximum is None or value <= maximum)):
            return float(value)
        raise ValueError(f"must be {text}")
    return check


def choice(values, null=None, ignore_case=False, null_value=None):
    """One of `values` (canonical spelling returned); null if `null` says what it means."""
    text = _or_null(f"one of {list(values)}", null)
    lookup = {v.lower() if ignore_case else v: v for v in values}

    def check(value):
        if value is None and null:
            return null_value
        key = value.lower() if (ignore_case and isinstance(value, str)) else value
        if isinstance(key, str) and key in lookup:
            return lookup[key]
        raise ValueError(f"must be {text}" + (" (any case)" if ignore_case else ""))
    return check


def text(null=None):
    """A string; null if `null` says what it means."""
    def check(value):
        if value is None and null:
            return None
        if isinstance(value, str) and value.strip():
            return value
        raise ValueError(f"must be {_or_null('a text', null)}")
    return check


def path():
    """A file path, or null."""
    return text(null="not used")


def paths():
    """A file path, a list of file paths, or null."""
    def check(value):
        if value is None:
            return None
        if isinstance(value, str) and value.strip():
            return value
        if isinstance(value, list) and value and all(isinstance(v, str) and v.strip() for v in value):
            return value
        raise ValueError("must be a file path, a list of file paths, or null")
    return check


def protein_forcefield(value):
    """amber, amber14, charmm or charmm_gromacs in any case; amber14-all and
    amber14-all.xml mean amber14, auto means amber (as the calculator reads them)."""
    names = {'amber': 'amber', 'amber14': 'amber14', 'charmm': 'charmm',
             'charmm_gromacs': 'charmm_gromacs', 'amber14-all': 'amber14',
             'amber14-all.xml': 'amber14', 'auto': 'amber'}
    if isinstance(value, str) and value.strip().lower() in names:
        return names[value.strip().lower()]
    raise ValueError("must be one of ['amber', 'amber14', 'charmm', 'charmm_gromacs'] (any case; "
                     "amber14-all and amber14-all.xml mean amber14)")


PLATFORMS = ('CUDA', 'OpenCL', 'CPU', 'Reference')


def canonical_platform(name):
    """Canonical OpenMM platform name of `name` in any case ('cuda' -> 'CUDA');
    None for null or 'auto' (automatic choice). Raises ValueError otherwise."""
    if name is None or (isinstance(name, str) and name.strip().lower() == 'auto'):
        return None
    if isinstance(name, str):
        for platform in PLATFORMS:
            if name.strip().lower() == platform.lower():
                return platform
    raise ValueError(f"must be one of {list(PLATFORMS)} (any case), or null/auto for the "
                     f"automatic choice (CUDA, then OpenCL, then CPU)")


def platform_name():
    return canonical_platform


# --------------------------------------------------------------------------
# The option table
# --------------------------------------------------------------------------

class _Required:
    def __repr__(self):
        return 'REQUIRED'


REQUIRED = _Required()


@dataclass(frozen=True)
class Option:
    """One configuration option.

    default  value used when the option is left out (REQUIRED: no default)
    check    value -> normalised value; raises ValueError saying what is allowed
    comment  its text in the generated template
    example  value shown in the template when it differs from the default
    heading  sub-heading written before it in the template
    path     True for input files that must exist
    """
    default: Any = None
    check: Optional[Callable[[Any], Any]] = None
    comment: str = ''
    example: Any = None
    heading: str = ''
    path: bool = False


ENTROPY_METHODS = ['none', 'interaction', 'c2', 'quasiharmonic', 'normal_mode', 'nmode']

OPTIONS: Dict[str, Dict[str, Option]] = {
    'input_files': {
        'complex_pdb': Option(
            REQUIRED, path(), path=True, example='path/to/complex.prmtop', comment=(
                "Native Mode (the topology of the MD simulation is read as is):\n"
                "  Amber .prmtop/.parm7, GROMACS .top (with a .gro of the same folder) or\n"
                "  .tpr, CHARMM/NAMD .psf (needs forcefield_settings.charmm_params), or a\n"
                "  serialised OpenMM System .xml with a .pdb of the same name next to it.\n"
                "Coordinate Mode (OpenGBSA builds the topology from a force field):\n"
                "  .pdb or .gro of the complex, plus ligand_mol for a small-molecule ligand.")),
        'trajectory': Option(
            REQUIRED, path(), path=True, example='path/to/trajectory.xtc',
            comment="MD trajectory (.xtc, .trr, .dcd, .nc, .trj, multi-model .pdb, ...)."),
        'solvated_topology': Option(None, path(), path=True, comment=(
            "Topology of the atoms in the trajectory when it still contains water and\n"
            "ions (for example the solvated .gro or .prmtop). Not needed for a dry\n"
            "trajectory.")),
        'receptor_topology': Option(None, path(), path=True, comment=(
            "Native Mode, Amber only (.prmtop or .parm7): receptor and ligand topologies\n"
            "written with the complex. Their atom counts define the receptor and the ligand.")),
        'ligand_topology': Option(None, path(), path=True),
        'ligand_mol': Option(None, path(), path=True, comment=(
            "Coordinate Mode: the small-molecule ligand (.sdf, or .mol2), with bond\n"
            "orders and, if possible, partial charges (else charge_method assigns them).\n"
            "Parameterised with OpenFF.")),
        'ligand_pdb': Option(None, path(), path=True, comment=(
            "Optional ligand-only .pdb; used only to check the ligand residue names.")),
        'ligand_resname': Option(None, text(null="detected"), comment=(
            "Residue name shown in the ligand box of the interactive 3D viewer.")),
    },
    'output_settings': {
        'output_directory': Option('mmgbsa_results', text(), comment=(
            "Each run writes <output_directory>/analysis_<date>_<time>/ (with _2, _3, ...\n"
            "added when that folder exists). opengbsa --output-dir DIR overrides this.")),
    },
    'analysis_settings': {
        # --- receptor and ligand ---
        'binding_mode': Option(
            'standard', choice(['standard', 'ppi']), heading='receptor and ligand', comment=(
                "standard: receptor = everything but the ligand, ligand found automatically\n"
                "  (receptor_topology atom count, else residue LIG/UNL/UNK/DRUG/MOL, else the\n"
                "  smallest non-solvent residue). ppi: receptor = chainid 0, ligand = chainid 1\n"
                "  of a topology with exactly two chains.")),
        'receptor_selection': Option(None, text(null="everything but the ligand"), comment=(
            "MDTraj selections; they override binding_mode, e.g. \"chainid 1\". With only\n"
            "ligand_selection the receptor is every other atom. The two must not overlap;\n"
            "atoms in neither are left out of the calculation.")),
        'ligand_selection': Option(None, text(null="detected")),
        'charge_method': Option('am1bcc', choice(['am1bcc', 'mmff94', 'gasteiger'], ignore_case=True),
                                comment=(
            "Coordinate Mode: partial charges of a ligand_mol file that has none, and of\n"
            "residues the protein force field lacks: am1bcc (needs AmberTools' sqm),\n"
            "mmff94 or gasteiger. Charges in the ligand file are always kept.")),
        # --- solvation model ---
        'gb_model': Option('OBC2', choice(['OBC2', 'OBC1', 'HCT', 'GBn', 'GBn2']),
                           heading='solvation model', comment="OBC2, OBC1, HCT, GBn, GBn2"),
        'gb_radii_set': Option('auto', choice(['auto', 'mbondi', 'mbondi2', 'mbondi3', 'bondi'],
                                              null="auto", null_value='auto'),
                               comment="auto, mbondi, mbondi2, mbondi3, bondi"),
        'salt_concentration': Option(0.15, number(minimum=0.0, maximum=2.0),
                                     comment="M, Debye-Hueckel screening"),
        'temperature': Option(300.0, number(greater_than=0.0, maximum=1000.0), comment="K"),
        'solute_dielectric': Option(1.0, number(greater_than=0.0),
                                    comment="interior dielectric of the GB and Coulomb terms"),
        'solvent_dielectric': Option(None, number(greater_than=0.0, null="78.5, or MMPBSA.py's "
                                                                       "value with amber_compatible"),
                                     comment="null = 78.5 (MMPBSA.py's value with amber_compatible)"),
        'sa_model': Option('LCPO', choice(['LCPO', 'ACE'], null="LCPO", null_value='LCPO'),
                           comment="non-polar term: LCPO or ACE"),
        'lcpo_typing': Option(None, choice(['element', 'amber'],
                                           null="element; amber with amber_compatible"),
                              comment="element (default) or amber"),
        'rgbmax': Option(None, number(greater_than=0.0, null="no Born-radius cutoff"),
                         comment="Born-radius cutoff in A; null = every atom pair"),
        'amber_compatible': Option(False, boolean(), comment=(
            "true = Amber MMPBSA.py's conventions: rgbmax 25 A, MMPBSA.py's solvent\n"
            "dielectric and LCPO with Amber's atom typing.")),
        # --- frames ---
        'frame_start': Option(None, integer(minimum=0, null="first frame"), heading='frames',
                              comment="0-based; null = first frame"),
        'frame_end': Option(None, integer(minimum=1, null="last frame"),
                            comment="exclusive; null = last frame"),
        'frame_stride': Option(None, integer(minimum=1, null="1"),
                               comment="every Nth frame of the range; null = 1"),
        'max_frames': Option(50, integer(minimum=1, maximum=10000, null="every frame"),
                             comment="null = every frame of the range after striding"),
        'frame_selection': Option('sequential', choice(['sequential', 'equidistant', 'random']), comment=(
            "Chooses max_frames frames among the frames of the range after striding:\n"
            "sequential: the first ones; equidistant: spread evenly; random: a random\n"
            "set (random_seed).")),
        'random_seed': Option(42, integer(null="42")),
        'reimage_trajectory': Option(True, boolean(),
                                     comment="make molecules whole across periodic boundaries"),
        # --- per-residue decomposition ---
        'run_per_residue_decomposition': Option(False, boolean(), heading='per-residue decomposition'),
        'decomp_frames': Option(10, integer(minimum=1, null="every analysed frame"), comment=(
            "frames decomposed: a subset of the analysed frames, chosen like\n"
            "frame_selection (the first ones, evenly spread, or random); not more than\n"
            "max_frames; null = every analysed frame")),
        'report_raw_energies': Option(False, boolean(),
                                      comment="also write complex/receptor/ligand energies"),
        'plot_top_residues': Option(10, integer(minimum=1),
                                    comment="residues in the hot-spot plots, table and report"),
        'energy_decomposition': Option(False, boolean(), comment=(
            "A second, faster decomposition inside the main calculation (fast or full),\n"
            "mainly for development.")),
        'decomposition_method': Option('full', choice(['full', 'fast'])),
        # --- entropy ---
        'run_entropy_analysis': Option(False, boolean(), heading='entropy (-TdS)', comment=(
            "run_entropy_analysis switches the estimate on; entropy_method chooses it:\n"
            "interaction (interaction entropy) or c2 (second-order cumulant), both from\n"
            "the fluctuation of the gas-phase interaction energy and reported with its\n"
            "standard deviation (IE unreliable above 3.6, C2 above 6.0 kcal/mol);\n"
            "quasiharmonic: ligand quasi-harmonic entropy (qha_analyze_complex: true\n"
            "adds receptor and complex to give a binding entropy); nmode: normal-mode\n"
            "binding entropy of complex, receptor and ligand, each minimized in GB\n"
            "without the non-polar term (amber_compatible: MMPBSA.py's HCT model,\n"
            "15 A Born-radius cutoff and no salt), on nmode_frames frames spread\n"
            "evenly over the analysed frames.")),
        'entropy_method': Option('interaction', choice(ENTROPY_METHODS)),
        'qha_analyze_complex': Option(False, boolean()),
        'nmode_frames': Option(10, integer(minimum=1, null="every analysed frame"),
                               comment="nmode: frames used; null = every analysed frame"),
        'nmode_drms': Option(0.001, number(greater_than=0.0),
                             comment="nmode: RMS gradient (kcal/mol/A) the minimization must reach"),
        'nmode_maxcyc': Option(10000, integer(minimum=1),
                               comment="nmode: iterations per minimization round"),
        'nmode_precision': Option('mixed', choice(['mixed', 'double'], null="double",
                                                  null_value='double'),
                                  comment="nmode: GPU force precision for the Hessian; double on CPU"),
        # --- execution ---
        'max_workers': Option(None, integer(minimum=1, null="min(cores - 1, 4)"), heading='execution',
                              comment="CPU worker processes; null = min(cores - 1, 4)"),
        'print_interval': Option(10, integer(minimum=1), comment="energy line every N frames"),
        # the former levels 2 and 3 behaved like 1 and are read as true
        'verbose': Option(True, boolean(allow_int=True), comment="false = fewer setup messages"),
    },
    'forcefield_settings': {
        'protein_forcefield': Option('amber', protein_forcefield, comment=(
            "Coordinate Mode only: amber (ff14SB, default), amber14 (amber14-all),\n"
            "charmm or charmm_gromacs (CHARMM36).")),
        'charmm_params': Option(None, paths(), comment=(
            "CHARMM/NAMD .psf input: parameter files (.prm, .str, .rtf), one path or a list.")),
        'charmm_coordinates': Option(None, path(), comment=(
            "CHARMM/NAMD .psf input: coordinates (.pdb, .coor, .crd) if no\n"
            "solvated_topology is given.")),
    },
    'platform_settings': {
        'preferred_platform': Option(None, canonical_platform, comment=(
            "null = automatic (CUDA, then OpenCL, then CPU). Or CUDA, OpenCL, CPU,\n"
            "Reference (any case). Select a GPU with the environment variable\n"
            "CUDA_VISIBLE_DEVICES.")),
        'decomposition_platform': Option(None, canonical_platform, comment=(
            "Platform for the per-residue decomposition; null = same as above.")),
        'strict_platform': Option(False, boolean(), comment=(
            "true = stop with an error instead of falling back to another platform, and\n"
            "stop when the per-residue decomposition fails.")),
    },
    'advanced_settings': {
        'n_jobs': Option(-1, integer(minimum=1, also=(-1,)), comment=(
            "Worker processes for a per-residue decomposition on the CPU; -1 = all cores.")),
    },
    'reporting_settings': {
        'include_all_plots': Option(True, boolean(), comment="energy time series and distributions"),
        'export_format': Option('png', choice(['png', 'svg', 'jpeg', 'webp']),
                                comment="format of the plot export button"),
        'export_width': Option(4200, integer(minimum=1)),
        'export_height': Option(2800, integer(minimum=1)),
    },
}

SECTION_TITLES = {
    'input_files': 'input files',
    'output_settings': 'output',
    'analysis_settings': 'analysis',
    'forcefield_settings': 'force field',
    'platform_settings': 'hardware',
    'advanced_settings': 'advanced',
    'reporting_settings': 'HTML report',
}

REQUIRED_SECTIONS = ('input_files', 'analysis_settings')

# Old spellings that still work, with the option that replaces them.
DEPRECATED_ALIASES = {
    ('params', 'output_directory'): ('output_settings', 'output_directory'),
    ('analysis_settings', 'per_residue_decomposition'): ('analysis_settings', 'run_per_residue_decomposition'),
    ('advanced_settings', 'run_per_residue_decomposition'): ('analysis_settings', 'run_per_residue_decomposition'),
}
# The legacy top-level 'input' section and its keys
LEGACY_INPUT_KEYS = {
    'topology': 'complex_pdb',
    'trajectory': 'trajectory',
    'ligand_mol': 'ligand_mol',
    'ligand_pdb': 'ligand_pdb',
    'receptor_topology': 'receptor_topology',
    'ligand_topology': 'ligand_topology',
    'solvated_topology': 'solvated_topology',
}

_RETIRED_RUNNER = "it was read only by the retired complete runner and had no effect on opengbsa runs"
REMOVED_SECTIONS = {
    'reproducibility_settings': _RETIRED_RUNNER,
    'visualization_settings': ("it was never read; the energy_decomposition 'fast' PyMOL script "
                               "(-0.5 kcal/mol, 5 A) and heat maps are always written"),
}
REMOVED_OPTIONS = {
    'input_files': {
        'receptor_resname': "it was never read",
        'ligand_charge': "it was never read",
    },
    'output_settings': {
        'analysis_name': _RETIRED_RUNNER,
    },
    'analysis_settings': {
        'use_cache': "the per-run cache never reused anything",
        'parallel_processing': ("CPU runs always split the frames over worker processes; "
                                "max_workers sets how many"),
        'dimer_mode': ("run a second analysis with receptor_selection and ligand_selection for "
                       "the interface between two units"),
        'receptor_subunits': ("run a second analysis with receptor_selection and ligand_selection "
                              "for the interface between two units"),
        'compound_name': _RETIRED_RUNNER,
        'save_frame_by_frame_csv': _RETIRED_RUNNER,
        'frame_by_frame_csv_name': _RETIRED_RUNNER,
        'include_residue_summary': _RETIRED_RUNNER,
        'frame_output_components': _RETIRED_RUNNER,
        'frame_output_format': _RETIRED_RUNNER,
    },
    'reporting_settings': {
        'entropy_approximation': "the report shows the result of entropy_method",
        'generate_final_report': _RETIRED_RUNNER,
        'include_plots': _RETIRED_RUNNER,
        'include_statistics': _RETIRED_RUNNER,
        'include_convergence': _RETIRED_RUNNER,
    },
}

NONBONDED_CUTOFF_MESSAGE = (
    "analysis_settings.nonbonded_cutoff is no longer supported: MM/GBSA evaluates every atom "
    "pair without a cutoff, as Amber's MMPBSA.py does (cut = 999). To limit the effective "
    "Born-radius integral, set analysis_settings.rgbmax (A).")


def known_options():
    """{section: [option, ...]} of every option in OPTIONS (for documentation checks)."""
    return {section: list(options) for section, options in OPTIONS.items()}


def default_config():
    """Every section and option with its default (REQUIRED options are None)."""
    return {section: {name: (None if opt.default is REQUIRED else copy.deepcopy(opt.default))
                      for name, opt in options.items()}
            for section, options in OPTIONS.items()}


def apply_aliases(config):
    """Map deprecated spellings in `config` (in place) onto the options that replace them.

    Returns the deprecation warnings. The legacy top-level 'input' section is
    copied into input_files when there is no input_files section;
    params.output_directory takes precedence over output_settings.output_directory,
    as before; the per_residue_decomposition aliases switch the decomposition
    on if any of the spellings does.
    """
    warnings = []
    if not isinstance(config, dict):
        return warnings
    # A section left empty in YAML loads as null
    for section, value in list(config.items()):
        if value is None and (section in OPTIONS or section in ('params', 'input')
                              or section.endswith(('_files', '_settings', 'params'))):
            config[section] = {}

    if 'input' in config:
        legacy = config.pop('input') or {}
        if 'input_files' in config:
            warnings.append("The top-level 'input' section is deprecated and is ignored because "
                            "an input_files section is given; use input_files only.")
        else:
            warnings.append("The top-level 'input' section is deprecated; use input_files "
                            "(its 'topology' key is input_files.complex_pdb).")
            config['input_files'] = {LEGACY_INPUT_KEYS.get(k, k): v for k, v in legacy.items()}

    for (old_section, old_key), (new_section, new_key) in DEPRECATED_ALIASES.items():
        old = config.get(old_section)
        if not isinstance(old, dict) or old_key not in old:
            continue
        value = old.pop(old_key)
        warnings.append(f"{old_section}.{old_key} is deprecated; use {new_section}.{new_key}.")
        target = config.setdefault(new_section, {})
        if not isinstance(target, dict):
            continue
        if new_key == 'run_per_residue_decomposition':
            target[new_key] = bool(target.get(new_key)) or bool(value)
        else:
            target[new_key] = value          # params.output_directory wins, as before
        if old_section == 'params' and not old:
            del config['params']

    settings = config.get('analysis_settings')
    if isinstance(settings, dict) and settings.get('binding_mode') == 'dimer_ligand':
        settings['binding_mode'] = 'standard'
        warnings.append("binding_mode 'dimer_ligand' is deprecated; it is the same as 'standard'.")
    return warnings


def effective_config(config):
    """The configuration a run uses: aliases mapped, every option present.

    Values are normalised where they pass their check (e.g. 'cuda' -> 'CUDA',
    300 -> 300.0); a value that fails its check is kept as given (the Python
    interface does not validate). Unknown and removed options are left out.
    """
    cfg = copy.deepcopy(config or {})
    apply_aliases(cfg)
    out = default_config()
    for section, options in OPTIONS.items():
        given = cfg.get(section)
        if not isinstance(given, dict):
            continue
        for name, opt in options.items():
            if name not in given:
                continue
            value = given[name]
            if opt.check is not None:
                try:
                    value = opt.check(value)
                except ValueError:
                    pass
            out[section][name] = value
    settings = out['analysis_settings']
    if settings.get('random_seed') is None:
        settings['random_seed'] = 42
    return out


def _suggestion(name, candidates):
    match = difflib.get_close_matches(name, list(candidates), n=1, cutoff=0.6)
    return f"; did you mean '{match[0]}'?" if match else ""


def _yaml_scalar(value):
    if value is None:
        return 'null'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, float):
        return repr(value)
    return str(value)


def complete_config_text():
    """The YAML text of `opengbsa --create-complete-config`, built from OPTIONS."""
    lines = [
        "# OpenGBSA complete configuration (opengbsa --create-complete-config).",
        "# Every option the analysis reads, set to its default. Delete what you do not",
        "# need: an option that is left out takes the value shown here. Options that",
        "# are not listed here are reported as unknown when the run starts.",
        "# Reference: https://opengbsa.readthedocs.io/",
        "",
    ]
    for section, options in OPTIONS.items():
        lines.append(f"# {'-' * 64} {SECTION_TITLES[section]}")
        lines.append(f"{section}:")
        for name, opt in options.items():
            if opt.heading:
                lines.append(f"  # --- {opt.heading} ---")
            value = opt.example if opt.example is not None else opt.default
            entry = f"  {name}: {_yaml_scalar(None if value is REQUIRED else value)}"
            comment = opt.comment.strip('\n')
            if comment and ('\n' in comment or len(entry) + len(comment) > 86):
                lines.extend(f"  # {line}".rstrip() for line in comment.split('\n'))
                lines.append(entry)
            elif comment:
                lines.append(f"{entry:<30} # {comment}")
            else:
                lines.append(entry)
        lines.append("")
    return "\n".join(lines)


class ConfigManager:
    """
    Configuration manager for MM/GBSA analysis.

    Handles loading, validation, and management of configuration files.
    """

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration manager.

        Args:
            config_path: Path to configuration file
        """
        self.config_path = config_path
        self.config = {}
        self.validation_errors = []
        self.validation_warnings = []

        if config_path:
            self.load_config(config_path)

    def load_config(self, config_path: str) -> bool:
        """
        Load configuration from YAML file.

        Args:
            config_path: Path to configuration file

        Returns:
            True if successful, False otherwise
        """
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self.config = yaml.safe_load(f) or {}

            # Normalize configuration (handle legacy 'input' section)
            self._normalize_config()

            self.config_path = config_path
            logger.info(f"Configuration loaded from: {config_path}")
            return True

        except FileNotFoundError:
            logger.error(f"Configuration file not found: {config_path}")
            return False
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML file: {e}")
            return False
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            return False

    def _normalize_config(self):
        """
        Map deprecated spellings (legacy 'input' section, params.output_directory,
        the per_residue_decomposition aliases) onto the current options and turn
        a section left empty in YAML (which loads as null) into an empty section.
        The deprecation warnings are reported by validate_config.
        """
        self._alias_warnings = apply_aliases(self.config)

    def validate_config(self) -> bool:
        """
        Validate configuration.

        Unknown sections and options give a warning (with the closest known
        name), removed options a warning that says why they are gone, and a
        value outside what an option accepts an error that says what it
        accepts. A valid configuration is completed with the default of every
        option it leaves out (see OPTIONS), so get_config() returns what the
        run uses.

        Returns:
            True if valid, False otherwise
        """
        self.validation_errors = []
        self.validation_warnings = list(getattr(self, '_alias_warnings', []))
        self._alias_warnings = []
        if not isinstance(self.config, dict):
            self.validation_errors.append("The configuration must be a YAML mapping of sections")
            return False
        self.validation_warnings += apply_aliases(self.config)

        # Check required sections
        for section in REQUIRED_SECTIONS:
            if section not in self.config:
                self.validation_errors.append(f"Missing required section: {section}")
        if self.validation_errors:
            return False

        self._check_sections()
        self._validate_input_files()
        self._validate_options()
        self._validate_cross_fields()

        if self.validation_errors:
            return False
        self.config = effective_config(self.config)
        return True

    def _check_sections(self):
        """Warn about unknown and removed sections and options."""
        for section in list(self.config):
            if section in OPTIONS:
                values = self.config[section]
                if values is None:
                    self.config[section] = {}
                elif not isinstance(values, dict):
                    self.validation_errors.append(f"Section {section} must be a mapping of options, "
                                                  f"not {values!r}")
                continue
            if section == 'params' and isinstance(self.config[section], dict):
                # params.output_directory was mapped above; nothing else was ever read
                for name in self.config[section]:
                    elsewhere = [s for s, o in OPTIONS.items() if name in o]
                    where = f"; it belongs in {elsewhere[0]}" if elsewhere else ""
                    self.validation_warnings.append(f"params.{name} is ignored (params only ever "
                                                    f"supplied output_directory){where}")
                continue
            if section in REMOVED_SECTIONS:
                self.validation_warnings.append(f"Section '{section}' was removed and is ignored: "
                                                f"{REMOVED_SECTIONS[section]}.")
            else:
                self.validation_warnings.append(f"Unknown section '{section}' is ignored"
                                                f"{_suggestion(section, OPTIONS)}")
        for section, options in OPTIONS.items():
            values = self.config.get(section)
            if not isinstance(values, dict):
                continue
            for name in values:
                if name in options:
                    continue
                if section == 'analysis_settings' and name == 'nonbonded_cutoff':
                    continue                           # checked in _validate_options
                removed = REMOVED_OPTIONS.get(section, {})
                if name in removed:
                    self.validation_warnings.append(f"Option '{name}' in {section} was removed and "
                                                    f"is ignored: {removed[name]}.")
                    continue
                hint = _suggestion(name, options)
                if not hint:
                    elsewhere = [s for s, o in OPTIONS.items() if name in o]
                    if elsewhere:
                        hint = f"; it belongs in {elsewhere[0]}"
                self.validation_warnings.append(f"Unknown option '{name}' in {section} is "
                                                f"ignored{hint}")

    def _validate_input_files(self):
        """Validate input files section."""
        input_files = self.config.get('input_files') or {}
        if not isinstance(input_files, dict):
            return
        for key, opt in OPTIONS['input_files'].items():
            if opt.default is REQUIRED and key not in input_files:
                self.validation_errors.append(f"Missing required file: {key}")
                continue
            if not opt.path or input_files.get(key) is None:
                continue
            file_path = input_files[key]
            if not isinstance(file_path, str):
                self.validation_errors.append(f"Invalid file path: input_files.{key} = {file_path!r} "
                                              f"must be a path")
            elif not self._validate_file_path(file_path, key):
                self.validation_errors.append(f"Invalid file path: {file_path} (input_files.{key}: "
                                              f"not found or not readable)")

    def _validate_options(self):
        """Check every given option against OPTIONS and store its normalised value."""
        for section, options in OPTIONS.items():
            values = self.config.get(section)
            if not isinstance(values, dict):
                continue
            for name, opt in options.items():
                if name not in values or opt.check is None or (section == 'input_files' and opt.path):
                    continue
                try:
                    values[name] = opt.check(values[name])
                except ValueError as e:
                    self.validation_errors.append(f"Invalid {section}.{name} = {values[name]!r}: {e}")
        settings = self.config.get('analysis_settings') or {}
        if 'nonbonded_cutoff' in settings:
            cutoff = settings.pop('nonbonded_cutoff')
            if cutoff is None or (_is_number(cutoff) and cutoff >= 999):
                self.validation_warnings.append(NONBONDED_CUTOFF_MESSAGE + " The value given means "
                                                "no cutoff, so the option is ignored.")
            else:
                self.validation_errors.append(NONBONDED_CUTOFF_MESSAGE)
        if 'entropy_method' in settings and settings['entropy_method'] is None:
            # choice() already reported it; say how to switch the estimate off
            self.validation_errors.append("entropy_method: null is not allowed; use 'none' or "
                                          "run_entropy_analysis: false for no entropy estimate")

    def _validate_cross_fields(self):
        """Validate cross-field dependencies (on the values given, before defaults)."""
        settings = self.config.get('analysis_settings', {})

        # Frame range validation
        frame_start = settings.get('frame_start')
        frame_end = settings.get('frame_end')
        if _is_int(frame_start) and _is_int(frame_end) and frame_end <= frame_start:
            self.validation_errors.append(f"frame_end ({frame_end}) must be greater than "
                                          f"frame_start ({frame_start}); frame_end is exclusive")

        # Decomposition frames validation (only when a decomposition runs)
        max_frames = settings.get('max_frames')
        decomp_frames = settings.get('decomp_frames')
        if (settings.get('run_per_residue_decomposition') and _is_int(decomp_frames)
                and _is_int(max_frames) and decomp_frames > max_frames):
            self.validation_errors.append(f"decomp_frames ({decomp_frames}) cannot be greater than "
                                          f"max_frames ({max_frames}): the decomposed frames are a "
                                          f"subset of the analysed frames")

        # Coordinate Mode builds a ligand chosen by a selection with the protein force
        # field, so a small molecule from ligand_mol would get no parameters
        input_files = self.config.get('input_files') or {}
        complex_pdb = str(input_files.get('complex_pdb') or '').lower()
        if (complex_pdb.endswith(('.pdb', '.gro')) and input_files.get('ligand_mol')
                and settings.get('ligand_selection')):
            self.validation_errors.append(
                "ligand_mol and ligand_selection cannot be combined in Coordinate Mode: a selected "
                "ligand is built with the protein force field. Leave out ligand_selection; the "
                "small-molecule ligand is then found by its residue name (LIG, UNL, UNK, DRUG, MOL)")

    def _validate_file_path(self, file_path: str, file_key: str) -> bool:
        """Validate file path."""
        if not file_path:
            return False

        # Check if file exists
        if not os.path.exists(file_path):
            return False

        # Check if file is readable
        if not os.access(file_path, os.R_OK):
            return False

        # Check file size (max 1GB)
        file_size = os.path.getsize(file_path)
        if file_size > 1024 * 1024 * 1024:  # 1GB
            self.validation_warnings.append(f"Large file: {file_path} ({file_size / 1024**3:.1f}GB)")

        return True

    def get_config(self) -> Dict[str, Any]:
        """Get configuration dictionary."""
        return self.config.copy()

    def get_input_files(self) -> Dict[str, str]:
        """Get input files configuration."""
        return self.config.get('input_files', {}).copy()

    def get_analysis_settings(self) -> Dict[str, Any]:
        """Get analysis settings configuration."""
        return self.config.get('analysis_settings', {}).copy()

    def get_validation_errors(self) -> List[str]:
        """Get validation errors."""
        return self.validation_errors.copy()

    def get_validation_warnings(self) -> List[str]:
        """Get validation warnings."""
        return self.validation_warnings.copy()

    def create_default_config(self, output_path: str) -> bool:
        """
        Create a short, commented configuration file with the settings most
        runs need (`opengbsa --create-config`).

        Args:
            output_path: Path to save configuration file

        Returns:
            True if successful, False otherwise
        """
        default_config_yaml = """# OpenGBSA configuration (opengbsa --create-config).
# Run it with:  opengbsa mmgbsa_config.yaml
# Every option is listed with --create-complete-config and explained at
# https://opengbsa.readthedocs.io/

input_files:
  # Native Mode: the topology used for the MD simulation (.prmtop, .top, .tpr,
  # .psf or a serialised OpenMM system .xml with a .pdb of the same name).
  # Coordinate Mode: a .pdb of the complex; OpenGBSA then builds the topology.
  complex_pdb: path/to/complex.prmtop
  trajectory: path/to/trajectory.xtc
  # Coordinate Mode only: the small-molecule ligand with bond orders and charges.
  # ligand_mol: path/to/ligand.sdf
  # Native Mode, optional: separate receptor and ligand topologies.
  # receptor_topology: path/to/receptor.prmtop
  # ligand_topology: path/to/ligand.prmtop
  # Topology of the atoms in the trajectory, if it still contains solvent.
  # solvated_topology: path/to/complex_solvated.gro

output_settings:
  output_directory: mmgbsa_results

analysis_settings:
  gb_model: OBC2              # OBC2, OBC1, HCT, GBn or GBn2
  salt_concentration: 0.15    # M
  temperature: 300.0          # K
  max_frames: 50
  frame_start: null           # first frame (0-based); null = start
  frame_end: null             # end frame (exclusive); null = last
  frame_stride: null          # every Nth frame; null = 1
  run_per_residue_decomposition: false
  decomp_frames: 10
  run_entropy_analysis: false
  entropy_method: interaction # interaction or c2
  amber_compatible: false     # true = reproduce Amber MMPBSA.py's conventions
"""

        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(default_config_yaml)

            logger.info(f"Default configuration created: {output_path}")
            return True

        except Exception as e:
            logger.error(f"Error creating default configuration: {e}")
            return False

    def create_complete_config(self, output_path: str) -> bool:
        """
        Create a configuration file with every option the analysis reads, its
        default and its meaning (`opengbsa --create-complete-config`), built
        from OPTIONS.

        Args:
            output_path: Path to save configuration file

        Returns:
            True if successful, False otherwise
        """
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(complete_config_text())

            logger.info(f"Complete configuration created: {output_path}")
            return True

        except Exception as e:
            logger.error(f"Error creating complete configuration: {e}")
            return False
