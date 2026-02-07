
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

for i, config_path in enumerate(configs):
    test_name = os.path.basename(config_path).replace('.yaml', '')
    print(f"[{i+1}/{len(configs)}] Running {test_name}...", end='', flush=True)
    
    # Run Command
    cmd = [
        "conda", "run", "-n", "mmgbsa", 
        "python", "-m", "mmgbsa.cli", 
        config_path
    ]
    
    try:
        # Capture output
        process = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
            check=False # Don't raise exception, handle return code
        )
        
        if process.returncode == 0:
            print(" PASS")
            passed += 1
            results.append((test_name, "PASS", ""))
        else:
            print(" FAIL")
            failed += 1
            # Get last 10 lines of error
            error_lines = process.stderr.strip().split('\n')[-10:]
            error_msg = "\n".join(error_lines)
            results.append((test_name, "FAIL", error_msg))
            
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
