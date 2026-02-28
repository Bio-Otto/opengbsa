
import os
import glob
import subprocess
import time
import sys

CONFIG_DIR = 'test/configs/comprehensive'
RESULTS_DIR = 'test/results/comprehensive'

# Find all yaml files
configs = sorted(glob.glob(os.path.join(CONFIG_DIR, '*.yaml')))

if not configs:
    print(f"No configuration files found in {CONFIG_DIR}")
    sys.exit(1)

print(f"Found {len(configs)} test configurations.")
print("-" * 60)

passed = 0
failed = 0
results = []

start_time = time.time()

# Optionally allow running a subset of configs by name prefix or full filename
selected = None
if len(sys.argv) > 1:
    # accept a single argument which may be the path or a glob
    selected = sys.argv[1]

for i, config_path in enumerate(configs):
    if selected:
        if selected not in config_path:
            continue
    test_name = os.path.basename(config_path).replace('.yaml', '')
    print(f"[{i+1}/{len(configs)}] Running {test_name}...", end='', flush=True)
    
    # (Previously tpr_custom was skipped when TprParser missing; the code
    # now handles fallback behavior so we run it unconditionally.)
    
    # Always use the currently running interpreter.
    # This avoids shell-level conda plugin side effects from `source activate ...`.
    cmd = [sys.executable, "-m", "mmgbsa.cli", config_path]
    
    try:
        # Capture output
        process = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
            check=False # Don't raise exception, handle return code
        )

        # Check for warning about ineffective selection
        combined = process.stdout + process.stderr
        warning_detected = "returned the same atom set as the default" in combined

        if process.returncode == 0:
            # Print PASS; but indicate if a warning was expected or found
            if warning_detected:
                print(" PASS (selection warning logged)")
            else:
                print(" PASS")
            passed += 1
            results.append((test_name, "PASS", ""))
            # Only the dedicated warning regression test requires this diagnostic.
            expect_warning = (test_name == "test_selection_warning")
            if expect_warning and not warning_detected:
                print(" WARNING: No selection warning detected for config that should have triggered one.")
                failed += 1
                results.append((test_name, "WARN_MISSING", "Missing expected warning"))
        else:
            print(" FAIL")
            failed += 1
            # Get last 10 lines of error
            error_lines = process.stderr.strip().split('\n')[-10:]
            error_msg = "\n".join(error_lines)
            results.append((test_name, "FAIL", error_msg))


        # Stop early if too many failures
        if failed >= 3:
            print("\nEarly exit: 3 or more failures encountered, halting further tests.")
            break

    except Exception as e:
        print(f" ERROR: {str(e)}")
        failed += 1
        results.append((test_name, "ERROR", str(e)))

end_time = time.time()
duration = end_time - start_time

print("-" * 60)
print(f"Testing Complete in {duration:.2f} seconds.")
print(f"Passed: {passed}")
print(f"Failed: {failed}")
print("-" * 60)

if failed > 0:
    print("\nFailures Details:")
    for name, status, msg in results:
        if status != "PASS":
            print(f"\nCould not run {name}:")
            print(msg)
            print("-" * 40)
    sys.exit(1)
else:
    print("\nAll tests passed successfully!")
    sys.exit(0)
