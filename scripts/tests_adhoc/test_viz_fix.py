import os
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Make sure we import local
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 mmgbsa
sys.path.insert(0, os.getcwd())

from mmgbsa.reporting import HTMLReportGenerator
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


output_dir = "test_viz_output"
os.makedirs(output_dir, exist_ok=True)

# Create dummy HTML file to simulate PandaMap
panda_html = os.path.join(output_dir, "structure_3d.html")
with open(panda_html, "w") as f:
    f.write("<html><body>PandaMap Content</body></html>")

generator = HTMLReportGenerator(output_dir)
# Dummy data
analysis_results = {'hot_spots': None}
frame_data = []

print("Testing PandaMap Case...")
# Generate report using HTML path (PandaMap case)
report_path_panda = generator.generate_report(analysis_results, frame_data, 
                                             complex_pdb_path=panda_html)

if report_path_panda:
    with open(report_path_panda, "r") as f:
        content = f.read()
        if 'id="structure-viewer" style="display:none;"' in content:
            print("PASS: Structure viewer is hidden for PandaMap.")
        else:
            print("FAIL: Structure viewer is NOT hidden for PandaMap.")
            # Debug: print relevant section
            idx = content.find('id="structure-viewer"')
            if idx != -1:
                print(f"Snippet: {content[idx:idx+100]}")
else:
    print("FAIL: Report generation failed for PandaMap case.")


print("\nTesting PDB Case...")
# Create dummy PDB
pdb_file = "test_dummy.pdb"
with open(pdb_file, "w") as f:
    f.write("ATOM  1  N   ALA A   1      1.000   1.000   1.000  1.00  0.00           N")

report_path_pdb = generator.generate_report(analysis_results, frame_data, complex_pdb_path=pdb_file)

if report_path_pdb:
    with open(report_path_pdb, "r") as f:
        content = f.read()
        # Should contain standard viewer class and NOT display:none
        if 'id="structure-viewer" class="viewer-container"' in content:
            print("PASS: Structure viewer is visible for PDB.")
        else:
            print("FAIL: Structure viewer logic incorrect for PDB.")
            idx = content.find('id="structure-viewer"')
            if idx != -1:
                print(f"Snippet: {content[idx:idx+100]}")
else:
    print("FAIL: Report generation failed for PDB case.")

# Clean up
try:
    os.remove(pdb_file)
    # import shutil
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # shutil.rmtree(output_dir)
except:
    pass
