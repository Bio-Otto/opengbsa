
import sys
import numpy as np

# Simulate the environment where numpy.compat is missing (NumPy 2.x)
# Effectively, we want to check if our patch in mmgbsa/core.py works.
# So we import mmgbsa.core first.

print("Importing mmgbsa.core...")
try:
    import mmgbsa.core
    print("✅ mmgbsa.core imported successfully")
except Exception as e:
    print(f"❌ Failed to import mmgbsa.core: {e}")
    sys.exit(1)

# Now check if numpy.compat exists in sys.modules
if 'numpy.compat' in sys.modules:
    print("✅ numpy.compat found in sys.modules")
else:
    print("❌ numpy.compat NOT found in sys.modules")
    sys.exit(1)

# Now try to import parmed (which triggers the error if not patched)
print("Importing parmed...")
try:
    import parmed
    print("✅ parmed imported successfully")
except ImportError as e:
    print(f"❌ Failed to import parmed: {e}")
    sys.exit(1)
except Exception as e:
    print(f"❌ An error occurred importing parmed: {e}")
    sys.exit(1)

print("🎉 Patch verification successful!")
