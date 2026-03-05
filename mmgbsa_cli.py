import mmgbsa; print(f'DEBUG: mmgbsa loaded from {mmgbsa.__file__}'); #!/usr/bin/env python3
"""
MM/GBSA CLI - Simple command line interface
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from mmgbsa.cli import main
except ImportError as e:
    print(f"Error importing mmgbsa: {e}")
    print("Make sure you have installed the package: pip install -e .")
    sys.exit(1)

if __name__ == "__main__":
    sys.exit(main()) 