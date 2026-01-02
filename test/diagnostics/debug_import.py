
import sys
from pathlib import Path
import os

print(f"Current working directory: {os.getcwd()}")
print(f"Sys path: {sys.path}")

try:
    sys.path.append(os.getcwd())
    from mmgbsa.complete_runner import CompleteMMGBSARunner
    print("Successfully imported CompleteMMGBSARunner")
except ImportError as e:
    print(f"Failed to import CompleteMMGBSARunner: {e}")
    import traceback
    traceback.print_exc()
except Exception as e:
    print(f"An error occurred: {e}")
    import traceback
    traceback.print_exc()
