"""
Command Line Interface for MM/GBSA Analysis Package

This module provides the command line interface for the MM/GBSA package.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional

try:
    from . import __version__, __author__, __email__
    from .config import ConfigManager
    from .runner import MMGBSARunner
    from .complete_runner import CompleteMMGBSARunner
    from .utils import setup_logging, create_output_directory, check_dependencies
except ImportError:
    # Fallback for direct execution
    __version__ = "0.0.1"
    __author__ = "H. Ibrahim Özdemir"
    __email__ = "halil.ibrahim.oozdemir@gmail.com"
    
    # Import modules directly
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    from mmgbsa.config import ConfigManager
    from mmgbsa.runner import MMGBSARunner
    from mmgbsa.complete_runner import CompleteMMGBSARunner
    from mmgbsa.utils import setup_logging, create_output_directory, check_dependencies

def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="MM/GBSA Analysis Package - Molecular Mechanics/Generalized Born Surface Area Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create default configuration
  mmgbsa --create-config

  # Run analysis with configuration file
  mmgbsa config.yaml

  # Run complete analysis
  mmgbsa --complete config.yaml

  # Check dependencies
  mmgbsa --check-deps

  # Show version
  mmgbsa --version
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
        action='store_true',
        help='Create default configuration file'
    )
    
    parser.add_argument(
        '--create-complete-config',
        action='store_true',
        help='Create complete configuration file with all options'
    )
    
    parser.add_argument(
        '--complete',
        action='store_true',
        help='Use complete runner with advanced features'
    )
    
    parser.add_argument(
        '--check-deps',
        action='store_true',
        help='Check if all dependencies are available'
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
        default='mmgbsa_results',
        help='Output directory (default: mmgbsa_results)'
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

def create_config_command(complete: bool = False):
    """Create configuration command."""
    config_name = "complete_mmgbsa_config.yaml" if complete else "mmgbsa_config.yaml"
    
    if os.path.exists(config_name):
        response = input(f"{config_name} already exists. Overwrite? (y/N): ")
        if response.lower() != 'y':
            print("Configuration creation cancelled.")
            return 0
    
    config_manager = ConfigManager()
    success = False
    
    if complete:
        success = config_manager.create_complete_config(config_name)
    else:
        success = config_manager.create_default_config(config_name)
    
    if success:
        print(f"✓ Configuration file created: {config_name}")
        print(f"Edit the file and run: mmgbsa {config_name}")
        return 0
    else:
        print(f"✗ Error creating configuration file: {config_name}")
        return 1

def run_analysis_command(config_file: str, complete: bool = False, output_dir: str = "mmgbsa_results"):
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
    
    # Create output directory
    try:
        output_path = create_output_directory(output_dir)
    except Exception as e:
        print(f"✗ Error creating output directory: {e}")
        return 1
    
    # Run analysis
    try:
        if complete:
            runner = CompleteMMGBSARunner(config_manager.get_config(), output_path)
        else:
            runner = MMGBSARunner(config_manager.get_config(), output_path)
        
        success = runner.run_analysis()
        
        if success:
            print(f"✓ Analysis completed successfully!")
            print(f"Results saved in: {output_path}")
            return 0
        else:
            print("✗ Analysis failed.")
            return 1
            
    except Exception as e:
        print(f"✗ Error during analysis: {e}")
        return 1

def show_version():
    """Show version information."""
    print(f"MM/GBSA Analysis Package v{__version__}")
    print(f"Author: {__author__}")
    print(f"Email: {__email__}")
    print(f"License: MIT")
    print(f"GitHub: https://github.com/halilibrahimozdemir/mmgbsa")

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
    
    if args.create_config:
        return create_config_command(complete=False)
    
    if args.create_complete_config:
        return create_config_command(complete=True)
    
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