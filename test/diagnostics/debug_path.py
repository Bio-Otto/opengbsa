
import os
import shutil
import sys

print(f"PATH: {os.environ.get('PATH')}")
print(f"antechamber location: {shutil.which('antechamber')}")
print(f"sqm location: {shutil.which('sqm')}")
print(f"python location: {sys.executable}")

try:
    from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
    print("Trying to instantiate AmberToolsToolkitWrapper...")
    wrapper = AmberToolsToolkitWrapper()
    print("Success!")
except Exception as e:
    print(f"Failed: {e}")
