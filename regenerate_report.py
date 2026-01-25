#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import yaml  # Added import
from pathlib import Path

# Add current directory to path to find mmgbsa package
sys.path.append(os.getcwd())

from mmgbsa.reporting import HTMLReportGenerator

def regenerate_report(result_dir):
    result_path = Path(result_dir)
    if not result_path.exists():
        print(f"Error: Directory {result_dir} not found.")
        return

    print(f"Regenerating report for: {result_dir}")
    
    # Load Config (for Plot control)
    config = {}
    try:
        config_path = "config_parallel_test.yaml"
        if os.path.exists(config_path):
             with open(config_path) as f:
                 config = yaml.safe_load(f)
             print(f"Loaded config from {config_path}")
        else:
             print("Warning: config_parallel_test.yaml not found. Using defaults.")
    except Exception as e:
         print(f"Warning: Failed to load config: {e}")

    # 1. Load Frame Data
    frame_csv = result_path / "frame_by_frame_decomposition.csv"
    if not frame_csv.exists():
        print("Error: frame_by_frame_decomposition.csv not found.")
        return
    
    print("Loading frame data...")
    df_frame = pd.read_csv(frame_csv)
    # Convert to list of dicts
    frame_data = df_frame.to_dict('records')
    
    # 2. Load Summary / Hotspots
    # We need per-residue data for hotspots
    hotspot_csv = result_path / "per_residue_detailed.csv"
    
    if hotspot_csv.exists():
        print(f"Loading hotspot data from {hotspot_csv.name}...")
        df_hotspots = pd.read_csv(hotspot_csv)
        
        # Check columns. Usually 'Residue', 'Total'
        # Normalize columns if needed
        df_hotspots.columns = [c.lower() for c in df_hotspots.columns]
        
        if 'total' in df_hotspots.columns:
            # Create 'residue_id' if not present (sometimes it is 'residue')
            if 'residue' in df_hotspots.columns and 'residue_id' not in df_hotspots.columns:
                 df_hotspots['residue_id'] = df_hotspots['residue']
            
            hot_spots = df_hotspots.sort_values('total', ascending=True).head(20)
        else:
            print("Warning: 'Total' column not found in per-residue CSV.")
            hot_spots = pd.DataFrame() # Empty
            
        analysis_results = {
            'mean_contribution': df_hotspots['total'].sum() if 'total' in df_hotspots.columns else 0.0,
            'n_residues': len(df_hotspots),
            'hot_spots': hot_spots
        }
    else:
        print("Warning: per_residue_detailed.csv not found. 3D coloring will be limited.")
        analysis_results = {'hot_spots': None}

    # 3. Locate PDB
    pdb_file = result_path / "temp_from_traj.pdb"
    if not pdb_file.exists():
        # Fallback to any pdb
        pdbs = list(result_path.glob("*.pdb"))
        if pdbs:
            pdb_file = pdbs[0]
        else:
            pdb_file = None
            
    print(f"Using PDB: {pdb_file}")
    
    # 4. Generate Report
    generator = HTMLReportGenerator(str(result_path))
    
    # Try to find ligand resname from config if possible, or infer?
    # For now, let's look for LIG inside the frame data columns?
    # Or just pass None and let it highlight generic hetflags.
    ligand_resname = None
    # Load Results Summary
    summary_file = result_path / 'results_summary.yaml'
    analysis_results = {}
    if summary_file.exists():
        try:
            import yaml
            with open(summary_file, 'r') as f:
                data = yaml.safe_load(f)
                if isinstance(data, list) and len(data) > 0:
                     analysis_results = data[0] # Take first item if list
                elif isinstance(data, dict):
                     analysis_results = data
        except Exception as e:
            print(f"Warning: Could not load results_summary.yaml: {e}")
    
    # Load Hotspots separately if not in summary
    hot_spots_df = pd.DataFrame()
            
    # Try to infer ligand name from PDB??
    if pdb_file:
         standard_res = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HOH', 'WAT', 'SOL', 'CL', 'NA', 'K', 'MG', 'KCX', 'CSO', 'PTR', 'SEP', 'TPO'}
         
         potential_ligands = set()
         with open(pdb_file, 'r') as f:
             for line in f:
                 if line.startswith("HETATM") or line.startswith("ATOM  "):
                     # PDB format: Name is cols 17-20
                     res_name = line[17:20].strip()
                     if res_name not in standard_res:
                         potential_ligands.add(res_name)
                         
         if potential_ligands:
             # Prioritize common ligand codes
             priority_list = ['LIG', 'UNL', 'MOL', 'DRG', 'INH']
             start_lig = None
             for p in priority_list:
                 if p in potential_ligands:
                     start_lig = p
                     break
             
             if not start_lig:
                 # Filter out likely modified residues if others exist
                 filtered = [l for l in potential_ligands if l not in ['KCX', 'CSO', 'PTR', 'TPO', 'SEP']]
                 if filtered:
                    start_lig = filtered[0]
                 else:
                    start_lig = list(potential_ligands)[0]
             
             ligand_resname = start_lig
             print(f"Potential ligands: {potential_ligands}. Selected: {ligand_resname}")
    
    print(f"Inferred Ligand Name: {ligand_resname}")
    
    # 2. Generate PandaMap 3D Report
    pandamap_html = None
    if pdb_file and ligand_resname:
        try:
            print("Generating PandaMap 3D Visualization...")
            from pandamap import HybridProtLigMapper
            from pandamap.create_3d_view import create_pandamap_3d_viz
            
            # Patch PDB (ATOM -> HETATM)
            patched_pdb = os.path.join(result_dir, "temp_fixed_for_panda.pdb")
            with open(pdb_file, 'r') as f_in, open(patched_pdb, 'w') as f_out:
                for line in f_in:
                    if line.startswith("ATOM") and (f" {ligand_resname} " in line):
                        line = "HETATM" + line[6:]
                    f_out.write(line)
            
            # Run PandaMap
            mapper = HybridProtLigMapper(patched_pdb, ligand_resname=ligand_resname)
            mapper.run_analysis()
            
            pandamap_html = os.path.join(result_dir, "structure_3d.html")
            create_pandamap_3d_viz(mapper, output_file=pandamap_html)
            print(f"PandaMap generated: {pandamap_html}")
            
            # Post-Process for Thinner Sticks & Hover Labels
            try:
                with open(pandamap_html, 'r') as f:
                    html_content = f.read()
                
                # 1. Thinner Sticks & Global Viewer Export
                # PandaMap defaults to 0.2 usually. We reduce it globally.
                html_content = html_content.replace("radius: 0.2", "radius: 0.14")
                html_content = html_content.replace("radius:0.2", "radius:0.14")
                
                # CRITICAL: Expose viewer globally so our script can find it
                html_content = html_content.replace("viewer = $3Dmol.createViewer", "window.viewer = viewer = $3Dmol.createViewer")
                
                # 2. Add Hover Labels (Robust) & Top 5 Permanent Labels & Ligand Style
                
                # Identify Top 5 Residues (from per_residue_detailed.csv)
                top_5_calls = ""
                try:
                    stats_csv = os.path.join(result_dir, "per_residue_detailed.csv")
                    if os.path.exists(stats_csv):
                         df_stats = pd.read_csv(stats_csv)
                         # Normalize columns
                         df_stats.columns = [c.lower() for c in df_stats.columns]
                         
                         if 'total' in df_stats.columns:
                             # Sort ascending (negative is better)
                             top_5 = df_stats.sort_values('total', ascending=True).head(5)
                             
                             for _, row in top_5.iterrows():
                                 # Logic to get number
                                 res_num = row.get('residue_number')
                                 if pd.isna(res_num): continue
                                 res_num = int(res_num)
                                 # OFFSET CORRECTION: CSV is +1 relative to PDB
                                 pdb_res_num = res_num - 1
                                 
                                 res_name = row.get('residue_name', 'RES')
                                 val = row['total']
                                 # Label Text: Use Original CSV Numbering (Matches Report Tables)
                                 label_text = f"{res_name}{res_num} ({val:.1f})"
                                 
                                 # Positioning: Use PDB Numbering (Matches Structure)
                                 top_5_calls += f"addSmartLabel(v, '{label_text}', {pdb_res_num});\n"
                except Exception as e:
                    print(f"Error calculating Top 5 items: {e}")

                hover_script = f"""
    <script>
      // Smart Label Function: Tries CA, then falls back to Residue Center
      function addSmartLabel(viewer, text, resNum) {{
          var selCA = {{resi: resNum, atom: 'CA'}};
          var atoms = viewer.getModel().selectedAtoms(selCA);
          var pos = (atoms.length > 0) ? atoms[0] : {{resi: resNum}};
          
          viewer.addLabel(text, {{
             position: pos,
             fontColor: 'black',
             fontSize: 14,
             showBackground: false,
             inFront: true
          }});
      }}

      function setupScene() {{
          var v = window.viewer || viewer;
          if (v) {{
              console.log("Setting up Top 5 Labels and Style...");
              // 1. PyMOL-like Ligand Style (Green Carbon)
              v.setStyle({{resn: '{ligand_resname}'}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.3}}}});
              
              // 2. Permanent Labels for Top 5 Residues
              try {{
                  {top_5_calls}
              }} catch(e) {{
                  console.error("Error adding top 5 labels:", e);
              }}
              
              // 3. Hover Labels
              v.setHoverable({{}}, true, 
                  function(atom, viewer, event, container) {{
                      if (!atom.label) {{
                          // Display Resi = PDB Resi + 1 (Offset Adjustment)
                          var displayResi = parseInt(atom.resi) + 1;
                          atom.label = viewer.addLabel(atom.resn + " " + displayResi, {{
                              position: atom, 
                              backgroundColor: 'rgba(0,0,0,0.7)', 
                              fontColor: 'white', 
                              fontSize: 12,
                              showBackground: true,
                              inFront: true,
                              backgroundOpacity: 0.7
                          }});
                      }}
                  }},
                  function(atom, viewer) {{
                      if (atom.label) {{
                          viewer.removeLabel(atom.label);
                          delete atom.label;
                      }}
                  }}
              );
              v.render();
          }} else {{
              console.log("Viewer not ready, retrying...");
              setTimeout(setupScene, 500);
          }}
      }}
      
      $(document).ready(function() {{
          // Wait slightly longer to ensure main script has run and assigned viewer
          setTimeout(setupScene, 2000); 
      }});
    </script>
    </body>
                """
                # Clean up previous injections by stripping old hover scripts if present
                # Not easy to regex out, but appending should override if logic is sound.
                # Just ensure we replace </body>
                     
                if "</body>" in html_content:
                    html_content = html_content.replace("</body>", hover_script)
                else:
                    html_content += hover_script
                    
                with open(pandamap_html, 'w') as f:
                    f.write(html_content)
                print("Applied styling and labels.")
                
            except Exception as e:
                print(f"Post-processing failed: {e}")
            
        except Exception as e:
            print(f"PandaMap generation failed: {e}")
            import traceback
            traceback.print_exc()

    # 3. Generate HTML Report
    print("Generating Interactive HTML Report...")
    output_html = os.path.join(result_dir, "interactive_report.html")
    
    # We pass 'pandamap_html' as the first argument to _generate_advanced_3dmol 
    # But HTMLReportGenerator logic might need adjustment if arguments shifted.
    # In reporting.py we changed _generate_advanced_3dmol signature to:
    # (pandamap_html_path, *args)
    # The call site in generate_report passes (complex_pdb_path, ligand_resname, analysis_results)
    # So we need to ensure reporting.py CALLER logic (generate_report method) is updated too?
    # Wait, I only updated the helper method. I need to check how generate_report calls it.
    
    # Let's assume for now I should instantiate and call generator manually or rely on existing flow?
    # Existing flow calls:
    # script, controls = self._generate_advanced_3dmol(complex_pdb_path, ligand_resname, analysis_results)
    #
    # My new signature: _generate_advanced_3dmol(pandamap_html_path, *args)
    # So complex_pdb_path argument will be treated as pandamap_html_path.
    # So I just need to pass the pandamap_html path as 'complex_pdb_path' when calling generate_report!
    # HACKY but efficient.
    
    # Load Global Results (for accurate plots)
    global_csv = result_path / "fixed_enhanced_mmgbsa_results_obc2.csv"
    global_results = None
    if global_csv.exists():
        try:
            df_g = pd.read_csv(global_csv)
            global_results = df_g.to_dict('records')
            print(f"Loaded Global Results: {len(global_results)} frames")
        except Exception as e:
            print(f"Failed to load global CSV: {e}")

    generator = HTMLReportGenerator(result_dir, config=config)
    # Correct order: analysis_results, frame_data, global_results
    report_path = generator.generate_report(analysis_results, frame_data, 
                              global_results=global_results,
                              complex_pdb_path=pandamap_html if pandamap_html else pdb_file, 
                              ligand_resname=ligand_resname)
    
    if report_path:
        print(f"\nSUCCESS! Report generated: {output_html}")
    else:
        print("\nFailed to generate report.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Auto-detect latest result directory
        dirs = [d for d in Path("mmgbsa_results").iterdir() if d.is_dir()]
        if dirs:
            latest = max(dirs, key=os.path.getmtime)
            print(f"No directory specified. Using latest: {latest}")
            regenerate_report(latest)
        else:
            print("Usage: python regenerate_report.py <result_directory>")
    else:
        regenerate_report(sys.argv[1])
