
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, AmberToolsToolkitWrapper, RDKitToolkitWrapper, OpenEyeToolkitWrapper

print("Registered toolkits:")
for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    print(f"  - {toolkit.__class__.__name__}: {toolkit.is_available()}")

print("\nSpecific checks:")
print(f"AmberToolsToolkitWrapper available: {AmberToolsToolkitWrapper().is_available()}")
print(f"RDKitToolkitWrapper available: {RDKitToolkitWrapper().is_available()}")
print(f"OpenEyeToolkitWrapper available: {OpenEyeToolkitWrapper().is_available()}")

import shutil
print(f"\nantechamber path: {shutil.which('antechamber')}")
