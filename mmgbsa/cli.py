"""
Command Line Interface for OpenGBSA

This module provides the command line interface for the OpenGBSA package.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional

from ._version import __version__, __author__, __email__, __link__, __license__

try:
    from .config import ConfigManager
    from .runner import MMGBSARunner, make_run_directory
    from .utils import setup_logging, check_dependencies
except ImportError:
    # Fallback for direct execution
    pass

DEFAULT_CONFIG_NAME = "mmgbsa_config.yaml"
COMPLETE_CONFIG_NAME = "complete_mmgbsa_config.yaml"


def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="OpenGBSA - Molecular Mechanics/Generalized Born Surface Area Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Write a configuration template: minimal, or with every option explained
  opengbsa --create-config
  opengbsa --create-complete-config my_config.yaml

  # Run an analysis defined by a configuration file
  opengbsa config.yaml
  opengbsa config.yaml --output-dir results

  # Check the installation
  opengbsa --version
  opengbsa --check-deps

  # List and download the test and validation data
  opengbsa --list-test-data
  opengbsa --fetch-test-data --only core-fixtures

Documentation: https://opengbsa.readthedocs.io/
        """
    )
    
    # Main arguments
    parser.add_argument(
        'config_file',
        nargs='?',
        help='Configuration file path (YAML format)'
    )
    
    # Options
    parser.add_argument(
        '--create-config',
        nargs='?', const=DEFAULT_CONFIG_NAME, default=None, metavar='PATH',
        help=f'Create a short configuration file (default PATH: {DEFAULT_CONFIG_NAME})'
    )
    
    parser.add_argument(
        '--create-complete-config',
        nargs='?', const=COMPLETE_CONFIG_NAME, default=None, metavar='PATH',
        help=f'Create a configuration file with every option and its default '
             f'(default PATH: {COMPLETE_CONFIG_NAME})'
    )
    
    parser.add_argument(
        '--complete',
        action='store_true',
        help='Deprecated and ignored: every run uses the same validated pipeline'
    )
    
    parser.add_argument(
        '--check-deps',
        action='store_true',
        help='Check if all dependencies are available'
    )

    parser.add_argument(
        '--fetch-test-data',
        action='store_true',
        help='Download and verify test/validation datasets from Zenodo (see mmgbsa/test_data_manifest.py)'
    )

    parser.add_argument(
        '--only',
        type=str,
        default=None,
        help='Comma-separated dataset names to fetch with --fetch-test-data (default: all)'
    )

    parser.add_argument(
        '--force-refetch',
        action='store_true',
        help='Re-download test data even if a valid cached copy exists (use with --fetch-test-data)'
    )

    parser.add_argument(
        '--list-test-data',
        action='store_true',
        help='List available test/validation datasets without downloading them'
    )
    
    parser.add_argument(
        '--version',
        action='store_true',
        help='Show version information'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--log-file',
        type=str,
        help='Log file path'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory (default: defined in config or mmgbsa_results)'
    )
    
    return parser

def check_dependencies_command():
    """Check dependencies command."""
    print("Checking MM/GBSA dependencies...")
    deps = check_dependencies()
    
    all_available = True
    for dep, available in deps.items():
        status = "✓" if available else "✗"
        print(f"  {status} {dep}")
        if not available:
            all_available = False
    
    if all_available:
        print("\n✓ All dependencies are available!")
        return 0
    else:
        print("\n✗ Some dependencies are missing.")
        print("Install missing dependencies with: pip install -r requirements.txt")
        return 1

def list_test_data_command():
    """List available test/validation datasets command."""
    from .fetch_test_data import list_datasets

    print("Available test/validation datasets:\n")
    for name, description, n_files, ready, size_gb, downloaded in list_datasets():
        if not ready:
            status = "not yet published (no checksums recorded)"
        else:
            status = "downloaded" if downloaded else "not downloaded"
        print(f"  {name} [{status}]")
        print(f"    {description}")
        size = f", {size_gb:.1f} GB" if size_gb else ""
        print(f"    {n_files} file(s){size}")
        print()
    print("Fetch with: opengbsa --fetch-test-data [--only name1,name2]")
    return 0

def fetch_test_data_command(only: Optional[str] = None, force: bool = False):
    """Download and verify test/validation datasets command."""
    from .fetch_test_data import data_root, fetch_all

    only_list = [n.strip() for n in only.split(',')] if only else None
    print(f"Fetching test data{f' ({only_list})' if only_list else ' (all datasets)'}...")
    try:
        results = fetch_all(only=only_list, force=force)
    except Exception as e:
        print(f"✗ Error fetching test data: {e}")
        return 1

    total_files = sum(len(files) for files in results.values())
    print(f"\n✓ Fetched {len(results)} dataset(s), {total_files} file(s) total, "
          f"available under {data_root() / 'test' / 'data'}")
    return 0

def create_config_command(complete: bool = False, path: Optional[str] = None):
    """Create configuration command: write the template to `path`."""
    config_name = path or (COMPLETE_CONFIG_NAME if complete else DEFAULT_CONFIG_NAME)
    
    if os.path.exists(config_name):
        try:
            response = input(f"{config_name} already exists. Overwrite? (y/N): ")
        except EOFError:
            response = ''
        if response.lower() != 'y':
            print("Configuration creation cancelled.")
            return 0
    parent = os.path.dirname(config_name)
    if parent:
        os.makedirs(parent, exist_ok=True)
    
    config_manager = ConfigManager()
    success = False
    
    if complete:
        success = config_manager.create_complete_config(config_name)
    else:
        success = config_manager.create_default_config(config_name)
    
    if success:
        print(f"✓ Configuration file created: {config_name}")
        print(f"Edit the file and run: opengbsa {config_name}")
        return 0
    else:
        print(f"✗ Error creating configuration file: {config_name}")
        return 1

def run_analysis_command(config_file: str, complete: bool = False, output_dir: Optional[str] = None):
    """Run analysis command."""
    print(f"Running MM/GBSA analysis with config: {config_file}")
    
    # Load and validate configuration
    config_manager = ConfigManager(config_file)
    
    if not config_manager.validate_config():
        print("✗ Configuration validation failed:")
        for error in config_manager.get_validation_errors():
            print(f"  - {error}")
        return 1
    
    # Show warnings
    warnings = config_manager.get_validation_warnings()
    if warnings:
        print("⚠ Configuration warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    
    # Create output directory (only if explicitly provided via CLI)
    output_path = None
    if output_dir is not None:
        try:
            output_path = str(make_run_directory(output_dir, 'mmgbsa_analysis'))
        except Exception as e:
            print(f"✗ Error creating output directory: {e}")
            return 1
    
    # Run analysis
    try:
        if complete:
            # The former "complete runner" ignored many settings (amber_compatible,
            # sa_model, entropy_method, ...) and is retired; one pipeline serves all runs.
            print("⚠️  --complete is deprecated and has no effect; running the standard pipeline.")
        runner = MMGBSARunner(config_manager.get_config(), output_path)
        
        success = runner.run_analysis()
        
        if success:
            print(f"✓ Analysis completed successfully!")
            print(f"Results saved in: {runner.output_dir}")
            return 0
        else:
            print("✗ Analysis failed.")
            return 1
            
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        return 1

def show_version():
    """Show version information."""
    print(f"OpenGBSA v{__version__}")
    print(f"Author: {__author__}")
    print(f"Email: {__email__}")
    print(f"License: {__license__}")
    print(f"GitHub: {__link__}")

def main():
    """Main CLI function."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level, args.log_file)
    
    # Handle special commands
    if args.version:
        show_version()
        return 0
    
    if args.check_deps:
        return check_dependencies_command()

    if args.list_test_data:
        return list_test_data_command()

    if args.fetch_test_data:
        return fetch_test_data_command(only=args.only, force=args.force_refetch)

    if args.create_config:
        return create_config_command(complete=False, path=args.create_config)
    
    if args.create_complete_config:
        return create_config_command(complete=True, path=args.create_complete_config)
    
    # Handle main analysis command
    if not args.config_file:
        parser.print_help()
        return 1
    
    if not os.path.exists(args.config_file):
        print(f"✗ Configuration file not found: {args.config_file}")
        return 1
    
    return run_analysis_command(
        args.config_file, 
        complete=args.complete,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    sys.exit(main()) 