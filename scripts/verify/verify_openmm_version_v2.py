import openmm
import openmm.app as app
import sys

print(f"OpenMM Version: {openmm.Platform.getOpenMMVersion()}")

try:
    if hasattr(openmm, 'LCPOForce'):
        print("openmm.LCPOForce is AVAILABLE")
    elif hasattr(app, 'LCPOForce'):
        print("openmm.app.LCPOForce is AVAILABLE")
    elif hasattr(app.internal, 'LCPOForce'):
        print("openmm.app.internal.LCPOForce is AVAILABLE")
    else:
        print("LCPOForce is NOT found")
except Exception as e:
    print(f"Error checking LCPO: {e}")
