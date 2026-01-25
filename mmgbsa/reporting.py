import os
import json
import logging
import pandas as pd
import numpy as np

# Configure logger
logger = logging.getLogger(__name__)

class HTMLReportGenerator:
    """
    Generates an interactive HTML report for MM/GBSA analysis results.
    Includes:
    - Summary statistics
    - Interactive Plotly Heatmap
    - Embedded 3D Structure Viewer (3Dmol.js)
    """
    
    def __init__(self, output_dir, config=None):
        self.output_dir = output_dir
        self.config = config or {}
        self.html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MM/GBSA Analysis Report</title>
    
    <!-- Bootstrap CSS -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    
    <!-- 3Dmol.js -->
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #f8f9fa; }
        .header-section { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 2rem 0; margin-bottom: 2rem; }
        .card { margin-bottom: 20px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border: none; }
        .card-header { background-color: #fff; border-bottom: 1px solid #eee; font-weight: bold; color: #333; }
        .viewer-container { height: 500px; width: 100%; position: relative; border: 1px solid #ddd; background: #fff; }
        .plot-container { height: 600px; width: 100%; }
        .stat-box { text-align: center; padding: 15px; }
        .stat-value { font-size: 24px; font-weight: bold; color: #2c3e50; }
        .stat-label { font-size: 14px; color: #7f8c8d; text-transform: uppercase; letter-spacing: 1px; }
    </style>
</head>
<body>

    <div class="header-section text-center">
        <div class="container">
            <h1>MM/GBSA Analysis Report</h1>
            <p class="lead">Interactive Detailed Analysis</p>
        </div>
    </div>

    <div class="container">
        <!-- Summary Statistics -->
        <div class="row mb-4">
            <div class="col-md-12">
                <div class="card">
                    <div class="card-header">Global Binding Statistics</div>
                    <div class="card-body">
                        <div class="row">
                            <!-- STATS_PLACEHOLDER -->
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- 3D Structure Viewer -->
        <div class="row mb-4">
            <div class="col-md-12">
                <div class="card">
                    <div class="card-header">
                        3D Structure Visualization
                        <span class="float-end text-muted small">Shift+Click to select • Drag to rotate</span>
                    </div>
                    <div class="card-body p-0">
                        <div id="structure-viewer" class="viewer-container"></div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Interactive Heatmap -->
        <div class="row mb-4">
            <div class="col-md-12">
                <div class="card">
                    <div class="card-header">Time-Evolution Interaction Heatmap</div>
                    <div class="card-body">
                         <div id="heatmap-plot" class="plot-container"></div>
                    </div>
                </div>
            </div>
        </div>
        
    </div>

    <div id="error-console" style="color:red; background:#ffe6e6; padding:10px; margin:20px; border:1px solid red; display:none;"></div>

    <!-- Initialization Scripts -->
    <script>
        function logError(msg) {
            var errDiv = document.getElementById('error-console');
            errDiv.style.display = 'block';
            errDiv.innerHTML += "<strong>Error:</strong> " + msg + "<br>";
            console.error(msg);
        }

        window.onerror = function(message, source, lineno, colno, error) {
            logError(message + " (" + source + ":" + lineno + ")");
        };

        try {
            // --- 3D Visualization Logic ---
            $(3D_VIZ_SCRIPT)
        } catch (e) {
            logError("3D Viewer Failed: " + e);
        }
        
        try {
            // --- Plotly Heatmap Logic ---
            $(PLOTLY_SCRIPT)
        } catch (e) {
            logError("Plotly Heatmap Failed: " + e);
        }
    </script>
</body>
</html>
"""

    
    def _calculate_entropy_term(self, df_glob, temp=300.0):
        """
        Calculate Entropic Contribution (-TdS) from Global Binding Energy fluctuations.
        Returns (neg_TdS_value, method_name)
        """
        try:
            reporting_settings = self.config.get('reporting_settings', {})
            method = reporting_settings.get('entropy_approximation', 'gaussian').lower()
            
            if method == 'none':
                return None, None
            
            if 'binding_energy' not in df_glob.columns:
                return None, None

            delta_H = df_glob['binding_energy'].values
            R = 0.001987 # kcal/mol/K
            RT = R * temp
            beta = 1.0 / RT
            
            # 1. Gaussian Approximation (Quasi-Harmonic)
            # Reliable for short simulations, assumes normal distribution of energies.
            # -TdS = sigma^2 / (2RT)
            sigma = np.std(delta_H, ddof=1) # Sample std deviation
            gaussian_val = (sigma ** 2) / (2 * RT)
            
            if method == 'gaussian':
                return gaussian_val, "Gaussian Approx"
            
            # 2. Interaction Entropy (Exponential Average)
            # -TdS = KT * ln < e^(beta * (E - <E>)) >
            # More rigorous but sensitive to outliers (convergence issues).
            elif method == 'exponential' or method == 'interaction_entropy':
                try:
                    delta_H_mean = np.mean(delta_H)
                    dev = delta_H - delta_H_mean
                    arg = beta * dev
                    
                    # LogSumExp for stability: log(mean(exp(x))) = log(sum(exp(x))) - log(N)
                    # We use numpy usually.
                    max_arg = np.max(arg)
                    # sum(exp(arg)) = exp(max_arg) * sum(exp(arg - max_arg))
                    # log(sum) = max_arg + log(sum(exp(arg - max_arg)))
                    sum_exp = np.sum(np.exp(arg - max_arg))
                    log_sum = max_arg + np.log(sum_exp)
                    log_mean = log_sum - np.log(len(arg))
                    
                    ie_val = RT * log_mean
                    return ie_val, "Interaction Entropy"
                except Exception as e:
                    logger.warning(f"Exponential entropy calculation failed: {e}")
                    return gaussian_val, "Gaussian Approx (Fallback)"
            
            return None, None
            
        except Exception as e:
            logger.error(f"Entropy calculation error: {e}")
            return None, None

    def generate_report(self, analysis_results, frame_data, global_results=None, complex_pdb_path=None, ligand_resname=None):
        """
        Generate the publication-ready interactive HTML report.
        """
        try:
            logger.info("Generating advanced HTML report...")
            import plotly.graph_objects as go
            import plotly.io as pio
            
            # 1. Prepare Base Template
            html_content = self.html_template.replace(
                '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>', 
                '<!-- Plotly Embedded -->'
            ).replace(
                '<script src="https://3Dmol.org/build/3Dmol-min.js"></script>',
                '<!-- 3Dmol Embedded -->'
            ).replace(
                '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">',
                '<style>/* Basic Reset */ body{margin:0;padding:20px;font-family:"Segoe UI",sans-serif;background:#f4f4f4;} .card{background:#fff;border-radius:8px;box-shadow:0 2px 10px rgba(0,0,0,0.1);margin-bottom:20px;padding:20px;} .row{display:flex;flex-wrap:wrap;gap:20px;} .col-md-12{width:100%;} .col-md-8{flex:2;} .col-md-4{flex:1;} .btn{padding:8px 16px;cursor:pointer;background:#007bff;color:#fff;border:none;border-radius:4px;} .btn:hover{background:#0056b3;} label{display:block;margin-top:10px;font-weight:600;} input[type=range]{width:100%;} .viewer_controls{background:#eee;padding:15px;border-radius:5px;} </style><!-- Bootstrap Removed -->'
            )
            
            # Error Console
            html_content = html_content.replace('<body>', '<body><div id="error-console" style="color:red; background:#ffe6e6; padding:10px; margin:20px; border:1px solid red; display:none;"></div>')
            
            # 2. Embed Libraries
            threedmol_js = self._get_3dmol_js()
            html_content = html_content.replace('<!-- 3Dmol Embedded -->', f'<script>{threedmol_js}</script>')
            
            # 3. Fill Statistics
            stats_html = self._generate_stats_html(analysis_results, global_results=global_results)
            html_content = html_content.replace("<!-- STATS_PLACEHOLDER -->", stats_html)
            
            # 4. Generate 3D Viewer Script (Advanced)
            viz_script, controls_html = self._generate_advanced_3dmol(complex_pdb_path, ligand_resname, analysis_results)
            html_content = html_content.replace("$(3D_VIZ_SCRIPT)", viz_script)
            
            # Inject PandaMap Iframe (full width)
            # We replace the placeholder div entirely with our content
            html_content = html_content.replace('<div id="structure-viewer" class="viewer-container"></div>', 
                                              f'<div style="width:100%;">{controls_html}</div>')
            
            # 5. Generate Plotly
            plotly_div = self._generate_plotly_div(frame_data, global_results=global_results)
            plotly_container = f'<div style="min-height: 700px; height: auto; width: 100%;">{plotly_div}</div>'
            html_content = html_content.replace('<div id="heatmap-plot" class="plot-container"></div>', plotly_container)
            html_content = html_content.replace('$(PLOTLY_SCRIPT)', '// Embedded Plotly')
            
            # Error Handler
            error_script = """<script>window.onerror = function(msg, url, line) { document.getElementById('error-console').style.display = 'block'; document.getElementById('error-console').innerHTML += 'Page Error: ' + msg + '<br>'; }</script>"""
            html_content = html_content.replace('</body>', f'{error_script}</body>')
            
            output_file = os.path.join(self.output_dir, "interactive_report.html")
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(html_content)
            return output_file
        except Exception as e:
            logger.error(f"Generate Report Failed: {e}")
            import traceback
            traceback.print_exc()
            return None
            
            # 6. Generate Plotly Heatmap (Embedded)
            # Include explicit height container to prevent collapse
            plotly_div = self._generate_plotly_div(frame_data)
            plotly_container = f'<div style="height: 700px; width: 100%;">{plotly_div}</div>'
            
            html_content = html_content.replace('<div id="heatmap-plot" class="plot-container"></div>', plotly_container)
            html_content = html_content.replace('$(PLOTLY_SCRIPT)', '// Embedded Plotly')
            
            # Add global error handler
            error_script = """
            <script>
            window.onerror = function(msg, url, line) {
                var c = document.getElementById('error-console');
                c.style.display = 'block';
                c.innerHTML += 'Error: ' + msg + ' (' + line + ')<br>';
            }
            </script>
            """
            html_content = html_content.replace('</body>', f'{error_script}</body>')
            
            # Save Report
            output_file = os.path.join(self.output_dir, "interactive_report.html")
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(html_content)
                
            logger.info(f"Report saved to: {output_file}")
            print(f"  Interactive report saved: {output_file}")
            return output_file
            
        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _get_3dmol_js(self):
        """Not used with PandaMap iframe."""
        return "// 3Dmol.js not needed (PandaMap embedded)"

    def _generate_plotly_div(self, frame_data, global_results=None):
        """Generate Multiple Plotly Divs (Dashboard) with Embedded JS."""
        if not frame_data:
            return "<div>No Data</div>"
        
        try:
            import pandas as pd
            import plotly.graph_objects as go
            import plotly.io as pio
            from plotly.subplots import make_subplots
            
            # Check Config
            reporting_settings = self.config.get('reporting_settings', {})
            include_all_plots = reporting_settings.get('include_all_plots', True)
            
            df = pd.DataFrame(frame_data)
            df_glob = pd.DataFrame(global_results) if global_results is not None else None
            
            # Data selection logic
            energy_cols = [col for col in df.columns if col.endswith('_total')]
            mean_energies = {col: df[col].mean() for col in energy_cols}
            top_n = 20
            sorted_residues = sorted(mean_energies.items(), key=lambda x: x[1])[:top_n]
            
            plots_html = []
            
            # --- Plot 1: Binding Energy Timeline (Components) ---
            if include_all_plots:
                # Component Logic (Hybrid approach if Global Data exists)
                
                # Decomp Sums (Approximation for VdW/Elec)
                vdw_cols = [c for c in df.columns if c.endswith('_vdw')]
                elec_cols = [c for c in df.columns if c.endswith('_electrostatic')]
                
                series_vdw = df[vdw_cols].sum(axis=1) if vdw_cols else pd.Series(0, index=df.index)
                series_elec = df[elec_cols].sum(axis=1) if elec_cols else pd.Series(0, index=df.index)
                
                # Global Data (Preferred for Total/Solvation due to self-energy/entropy)
                if df_glob is not None and not df_glob.empty:
                    series_total = df_glob['binding_energy']
                    series_solv = df_glob['delta_gb'] + df_glob['delta_sa']
                else:
                    # Fallback (Warn: Solvation incorrect)
                    solv_cols = [c for c in df.columns if c.endswith('_solvation')]
                    series_solv = df[solv_cols].sum(axis=1) if solv_cols else pd.Series(0, index=df.index)
                    series_total = df[energy_cols].sum(axis=1)
                
                # Line Plot
                # Line Plot
                fig_line = go.Figure()
                fig_line.add_trace(go.Scatter(x=df['frame_number'], y=series_total, mode='lines', name='Total Binding (ΔH)', line=dict(color='black', width=3)))
                
                # Add Free Energy (Evidence of Entropy)
                neg_TdS, s_method = self._calculate_entropy_term(df_glob)
                if neg_TdS is not None:
                     series_free = series_total + neg_TdS
                     fig_line.add_trace(go.Scatter(x=df['frame_number'], y=series_free, mode='lines', name='Est. Free Energy (ΔG)', line=dict(color='#d62728', width=2.5)))
                     
                fig_line.add_trace(go.Scatter(x=df['frame_number'], y=series_vdw, mode='lines', name='Van der Waals', line=dict(color='green', width=1.5, dash='dot')))
                fig_line.add_trace(go.Scatter(x=df['frame_number'], y=series_elec, mode='lines', name='Electrostatic', line=dict(color='blue', width=1.5, dash='dot')))
                fig_line.add_trace(go.Scatter(x=df['frame_number'], y=series_solv, mode='lines', name='Solvation', line=dict(color='orange', width=1.5, dash='dot')))

                fig_line.update_layout(title='Energy Components vs Time (inc. Entropy)', xaxis_title='Frame', yaxis_title='Energy (kcal/mol)', height=500, legend=dict(orientation="h", y=1.1))
                plots_html.append(pio.to_html(fig_line, full_html=False, include_plotlyjs=True)) 

                # --- Plot 2: Mean Components Bar Chart (Summary) ---
                means = {
                    'Van der Waals': series_vdw.mean(),
                    'Electrostatic': series_elec.mean(),
                    'Solvation': series_solv.mean(),
                    'Total Binding': series_total.mean()
                }
                sems = {
                    'Van der Waals': series_vdw.sem(),
                    'Electrostatic': series_elec.sem(),
                    'Solvation': series_solv.sem(),
                    'Total Binding': series_total.sem()
                }
                
                fig_summary = go.Figure(go.Bar(
                    x=list(means.keys()),
                    y=list(means.values()),
                    error_y=dict(type='data', array=list(sems.values())), # Add Error Bars (SEM)
                    marker_color=['green', 'blue', 'orange', 'black']
                ))
                fig_summary.update_layout(title='Average Energy Components (kcal/mol)', yaxis_title='Energy (kcal/mol)', height=500)
                plots_html.append(pio.to_html(fig_summary, full_html=False, include_plotlyjs=False))

                # --- Plot 3: Top Residues (Bar Chart) ---
                bar_names = [col.replace('_total', '').replace('_A', '') for col, _ in sorted_residues]
                bar_values = [val for _, val in sorted_residues]
                bar_names = bar_names[::-1]
                bar_values = bar_values[::-1]
                
                fig_bar = go.Figure(go.Bar(x=bar_values, y=bar_names, orientation='h', marker=dict(color=bar_values, colorscale='RdBu_r', cmin=-5, cmax=5)))
                fig_bar.update_layout(title=f'Top {top_n} Contributing Residues', xaxis_title='Mean Interaction Energy (kcal/mol)', height=500)
                plots_html.append(pio.to_html(fig_bar, full_html=False, include_plotlyjs=False))

            # --- Plot 4: Heatmap (Always Included) ---
            z_data = []
            y_labels = []
            x_labels = df['frame_number'].tolist()
            
            for col, mean_val in sorted_residues:
                clean_name = col.replace('_total', '').replace('_A', '')
                y_labels.append(f"{clean_name} ({mean_val:.1f})")
                z_data.append(df[col].tolist())
            
            fig_heat = go.Figure(data=go.Heatmap(
                z=z_data,
                x=x_labels,
                y=y_labels,
                colorscale='RdBu',
                reversescale=True,
                zmid=0,
                zmin=-5.0, 
                zmax=5.0,
                colorbar=dict(title='Energy')
            ))
            
            fig_heat.update_layout(
                title='Per-Residue Interaction Energy Timeline',
                xaxis_title='Frame Number',
                yaxis_title='Residue',
                autosize=True,
                height=600,
                margin=dict(l=150, r=50, b=50, t=50)
            )
            
            # If line/bar plots were skipped, this is the first plot, so include JS
            include_js = True if not plots_html else False
            plots_html.append(pio.to_html(fig_heat, full_html=False, include_plotlyjs=include_js))

            # --- Plot 5: Binding Energy Distribution (Histogram) ---
            if include_all_plots and df_glob is not None:
                fig_dist = go.Figure()
                
                # Trace 1: Binding Energy (Enthalpy)
                fig_dist.add_trace(go.Histogram(
                    x=series_total,
                    histnorm='probability density',
                    name='Binding Energy (ΔH)',
                    marker_color='#6f42c1',
                    opacity=0.6
                ))
                
                # Trace 2: Free Energy (Enthalpy + Entropy)
                if neg_TdS is not None:
                     series_free = series_total + neg_TdS
                     fig_dist.add_trace(go.Histogram(
                        x=series_free,
                        histnorm='probability density',
                        name='Est. Free Energy (ΔG)',
                        marker_color='#d62728',
                        opacity=0.6
                    ))
                     # Add Mean Lines for both
                     mean_G = series_free.mean()
                     fig_dist.add_vline(x=mean_G, line_width=2, line_dash="dash", line_color="#d62728")
                     fig_dist.add_annotation(x=mean_G, y=1.05, yref="paper", text=f"ΔG: {mean_G:.1f}", showarrow=False, bgcolor="white")

                # Add Mean Line for Enthalpy
                mean_H = series_total.mean()
                fig_dist.add_vline(x=mean_H, line_width=2, line_dash="dash", line_color="#6f42c1")
                fig_dist.add_annotation(x=mean_H, y=1.05, yref="paper", text=f"ΔH: {mean_H:.1f}", showarrow=False, bgcolor="white")

                fig_dist.update_layout(
                    title='Distribution: Binding Energy vs Free Energy',
                    xaxis_title='Energy (kcal/mol)',
                    yaxis_title='Density',
                    height=500,
                    bargap=0.1,
                    barmode='overlay' # Crucial for overlapping histograms
                )
                plots_html.append(pio.to_html(fig_dist, full_html=False, include_plotlyjs=False))

            # Concatenate
            dashboard_html = '<div class="row">'
            for plot in plots_html:
                dashboard_html += f'<div class="col-12 mb-4">{plot}</div>'
            dashboard_html += '</div>'
            
            return dashboard_html
            
        except Exception as e:
            return f"<div>Plot Generation Failed: {e}</div>"

    def _generate_stats_html(self, results):
        """Generate HTML for summary statistics cards."""
        if not results:
            return "<div class='alert alert-warning'>No results available</div>"
            
        mean_bind = results.get('total_contribution', 0.0) # Note: 'total_contribution' in analysis_results refers to specific sum, verify structure
        # Actually analysis_results structure from decomposition.py:
        # 'total_contribution' = df['total'].sum() for one frame? No, it's averaged results.
        # Let's use simple generic stats if available
        
        # We prefer Global Binding Energy if possible, but 'analysis_results' passed here is usually the Per-Residue summary.
        # Let's assume we pass the global 'mmgbsa_results' mainly? 
        # Or better, we just show Hotspot info.
        
    def _generate_stats_html(self, analysis_results, global_results=None):
        """Generate Global Statistics HTML (Table)."""
        html = ""
        
        # Determine Source Data
        mean_bind = "N/A"
        std_bind = "N/A"
        sem_bind = "N/A"
        
        df = None
        if global_results:
             import pandas as pd
             df = pd.DataFrame(global_results)
             
             if 'binding_energy' in df.columns:
                 mean_bind = f"{df['binding_energy'].mean():.2f}"
                 std_bind = f"{df['binding_energy'].std():.2f}"
                 sem_bind = f"{df['binding_energy'].sem():.2f}"
             
             if 'delta_nb' in df.columns:
                 mean_nb = df['delta_nb'].mean()
                 mean_gb = df['delta_gb'].mean()
                 mean_sa = df['delta_sa'].mean()
        
        # Create Enhanced Table
        html += """
        <div class="col-md-12">
            <table class="table table-bordered table-hover" style="background:white; border-radius:8px; overflow:hidden;">
                <thead class="table-light">
                    <tr>
                        <th>Metric</th>
                        <th>Mean (kcal/mol)</th>
                        <th>Std Dev</th>
                        <th>SEM</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Total Binding
        html += f"""
                    <tr style="font-weight:bold; background:#f8f9fa;">
                        <td>Total Binding Energy (ΔG)</td>
                        <td>{mean_bind}</td>
                        <td>{std_bind}</td>
                        <td>{sem_bind}</td>
                    </tr>
        """
        
        # Components (if DF available)
        if df is not None:
             # Calculate Entropy using helper
             neg_TdS, method_name = self._calculate_entropy_term(df)
             
             if neg_TdS is not None:
                 delta_H_mean = df['binding_energy'].mean()
                 delta_G_total = delta_H_mean + neg_TdS
                 
                 html += f"""
                    <tr>
                        <td>Entropic Contribution (-TΔS) <br><small class="text-muted">({method_name})</small></td>
                        <td>{neg_TdS:.2f}</td>
                        <td>N/A</td>
                        <td>N/A</td>
                    </tr>
                    <tr style="background:#e8f4f8; font-weight:bold;">
                        <td>Total ΔG (inc. Entropy)</td>
                        <td>{delta_G_total:.2f}</td>
                        <td>-</td>
                        <td>-</td>
                    </tr>
                 """
             elif self.config.get('reporting_settings', {}).get('entropy_approximation') == 'none':
                 pass # User explicitly disabled
             else:
                 pass # Approximation not available (e.g. no binding energy column)

             # Solvation
             s_mean = df['delta_gb'] + df['delta_sa']
             html += f"""
                    <tr>
                        <td>Solvation Energy (ΔG_solv)</td>
                        <td>{s_mean.mean():.2f}</td>
                        <td>{s_mean.std():.2f}</td>
                        <td>{s_mean.sem():.2f}</td>
                    </tr>
             """
             # Electrostatic / VdW (NB)
             if 'delta_nb' in df.columns:
                 nb = df['delta_nb']
                 html += f"""
                    <tr>
                        <td>Non-Bonded (VdW + Elec)</td>
                        <td>{nb.mean():.2f}</td>
                        <td>{nb.std():.2f}</td>
                        <td>{nb.sem():.2f}</td>
                    </tr>
                 """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        # If no global data, fallback to old box style for simple generic stats
        if not global_results and analysis_results:
             # Fallback code
             mean_contr = analysis_results.get('mean_contribution', 0)
             html += f"<div class='alert alert-secondary'>Detailed global statistics unavailable. Mean Residue Contribution: {mean_contr:.2f}</div>"
             
        return html

    
    def _generate_advanced_3dmol(self, pandamap_html_path, *args):
        """Embed PandaMap HTML output."""
        if not pandamap_html_path or not os.path.exists(pandamap_html_path):
            return "// No PandaMap", "<div class='alert alert-warning'>PandaMap 3D report not found.</div>"

        with open(pandamap_html_path, 'r', encoding='utf-8') as f:
            panda_content = f.read()
        
        # Escape for srcdoc
        import html
        escaped_content = html.escape(panda_content)
        
        controls_html = f"""
        <div class="alert alert-info py-1 mb-2">
            <small><b>PandaMap Integration:</b> Visualization generated by PandaMap package.</small>
        </div>
        <!-- Embedding full PandaMap HTML in iframe -->
        <iframe srcdoc="{escaped_content}" style="width: 100%; height: 750px; border: none; border-radius: 8px;"></iframe>
        """
        
        # No extra script needed as iframe handles it
        script = "// PandaMap embedded via iframe"
        
        return script, controls_html
        
        # Build Hotspot JSON
        hotspots_json = "[]"
        if analysis_results and 'hot_spots' in analysis_results:
            df = analysis_results['hot_spots']
            hotspots_list = []
            
            # Robust extraction logic from before
            if not df.empty and 'residue_number' in df.columns:
                for _, row in df.iterrows():
                    hotspots_list.append({
                        'resi': int(row['residue_number']),
                        'resn': row['residue_name'],
                        'chain': 'A',
                        'energy': float(row['total'])
                    })
            elif not df.empty and 'residue_id' in df.columns:
                 import re
                 for _, row in df.iterrows():
                     rid = str(row['residue_id'])
                     parts = rid.split('_')
                     resi = None
                     resn = ""
                     if len(parts) >= 2 and parts[1].isdigit():
                         resn = parts[0]
                         resi = int(parts[1])
                     else:
                         match = re.search(r'(\d+)', rid)
                         if match: resi = int(match.group(1))
                     
                     if resi:
                         hotspots_list.append({
                            'resi': resi,
                            'resn': resn, # Pass name too
                            'energy': float(row['total'])
                         })
            import json
            hotspots_json = json.dumps(hotspots_list)

        # Build Interaction Table HTML
        table_rows = ""
        if analysis_results and 'hot_spots' in analysis_results:
             df = analysis_results['hot_spots']
             # Reuse hotspots_list from JSON generation if possible, or iterate again
             # Let's iterate the 'hotspots_list' logic we built for JS json
             # We can't access it easily here unless we reused the list. 
             # Re-building list for python-side table generation:
             top_int = []
             if not df.empty and 'residue_number' in df.columns:
                 for _, row in df.iterrows():
                     top_int.append((row['residue_name'], row['residue_number'], row['total']))
             elif not df.empty and 'residue_id' in df.columns:
                 import re
                 for _, row in df.iterrows():
                     parts = str(row['residue_id']).split('_')
                     if len(parts)>=2 and parts[1].isdigit():
                         top_int.append((parts[0], parts[1], row['total']))
                     else:
                         match = re.search(r'(\d+)', str(row['residue_id']))
                         if match: top_int.append(("", match.group(1), row['total']))
            
             for resn, resi, ene in top_int:
                 res_label = f"{resn} {resi}"
                 color_class = "text-danger" if ene < -2.0 else "text-dark"
                 table_rows += f"<tr><td><span class='badge bg-light text-dark border'>{res_label}</span></td><td class='{color_class} fw-bold'>{ene:.2f}</td></tr>"

        # Controls HTML + Table
        controls_html = f"""
        <div class="mb-3">
            <label>Ligand Residue Name</label>
            <input type="text" id="ligandInput" class="form-control form-control-sm" value="{ligand_resname or ''}" onchange="updateViz()">
        </div>

        <ul class="nav nav-tabs" id="vizTabs">
          <li class="nav-item"><a class="nav-link active" data-bs-toggle="tab" href="#controls">Controls</a></li>
          <li class="nav-item"><a class="nav-link" data-bs-toggle="tab" href="#table">Interactions</a></li>
        </ul>

        <div class="tab-content mt-2">
            <div class="tab-pane fade show active" id="controls">
                <label>Style</label>
                <select id="styleSelect" onchange="updateViz()" class="form-select form-select-sm">
                    <option value="cartoon">Cartoon</option>
                    <option value="surface">Molecular Surface</option>
                    <option value="stick">Stick</option>
                </select>
                
                <label>Coloring</label>
                <select id="colorSelect" onchange="updateViz()" class="form-select form-select-sm">
                    <option value="energy">Binding Energy</option>
                    <option value="spectrum">Rainbow</option>
                    <option value="hydro">Hydrophobicity</option>
                </select>
                
                <label class="mt-2 text-muted text-uppercase" style="font-size:0.75rem; font-weight:bold;">Toggles</label>
                <div class="form-check form-switch">
                  <input class="form-check-input" type="checkbox" id="showHbonds" checked onchange="updateViz()">
                  <label class="form-check-label" for="showHbonds">H-Bonds (Struct)</label>
                </div>
                <div class="form-check form-switch">
                  <input class="form-check-input" type="checkbox" id="showEnergeticLinks" checked onchange="updateViz()">
                  <label class="form-check-label" for="showEnergeticLinks">Energy Links</label>
                </div>
                <div class="form-check form-switch">
                    <input class="form-check-input" type="checkbox" id="showLigand" checked onchange="updateViz()">
                    <label class="form-check-label" for="showLigand">Ligand</label>
                </div>
                <div class="form-check form-switch">
                    <input class="form-check-input" type="checkbox" id="showLabels" checked onchange="updateViz()">
                    <label class="form-check-label" for="showLabels">Labels</label>
                </div>

                <div class="d-grid gap-2 mt-3">
                    <button class="btn btn-primary btn-sm" onclick="takeSnapshot()">📷 Snapshot</button>
                    <button class="btn btn-outline-secondary btn-sm" onclick="viewer.zoomTo()">⟲ Reset</button>
                </div>
            </div>
            
            <div class="tab-pane fade" id="table">
                <div class="table-responsive" style="max_height: 300px; overflow-y: auto;">
                    <table class="table table-sm table-hover text-center" style="font-size:0.9rem;">
                        <thead class="table-light"><tr><th>Residue</th><th>Energy (kcal)</th></tr></thead>
                        <tbody>
                            {table_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        """
        
        # JS Logic
        script = f"""
        let element = document.getElementById('structure-viewer');
        let config = {{ backgroundColor: 'white', antialias: true }};
        window.viewer = $3Dmol.createViewer(element, config);
        
        // Add Bootstrap JS for Tabs if missing (minimal wrapper)
        // Actually tabs are CSS/JS dependent. Since we removed Bootstrap JS, we need simple JS for tabs?
        // Or we just show one by default.
        // Let's add a simple script to handle tab clicks since libraries are offline.
        
        let pdbData = `{pdb_data}`;
        let hotspots = {hotspots_json};
        
        window.viewer.addModel(pdbData, "pdb");
        
        window.updateViz = function() {{
            let v = window.viewer;
            v.removeAllShapes(); // Removes lines/cylinders
            v.removeAllSurfaces();
            v.removeAllLabels();
            v.setStyle({{}});
            
            let ligVal = document.getElementById('ligandInput').value.trim();
            let style = document.getElementById('styleSelect').value;
            let colorScheme = document.getElementById('colorSelect').value;
            let showHbonds = document.getElementById('showHbonds').checked;
            let showLinks = document.getElementById('showEnergeticLinks').checked;
            let showLigand = document.getElementById('showLigand').checked;
            let showLabels = document.getElementById('showLabels').checked;
            
            // Selections
            let proteinSel = ligVal ? {{not: {{resn: ligVal}}}} : {{hetflag: false}};
            proteinSel = {{and: [proteinSel, {{not: {{resn: ['HOH','WAT','SOL','CL','NA','MG']}} }}]}};
            let ligandSel = ligVal ? {{resn: ligVal}} : {{hetflag: true}};
            
            // 1. Protein Style
            let pStyle = {{}};
            if (colorScheme === 'spectrum') pStyle = {{color: 'spectrum'}};
            else if (colorScheme === 'hydro') pStyle = {{colorscheme: 'amino'}};
            else if (colorScheme === 'chain') pStyle = {{colorscheme: 'chain'}};
            else pStyle = {{color: 'white'}}; 
            
            if (style === 'cartoon') v.setStyle(proteinSel, {{cartoon: {{style:'oval', ...pStyle}}, stick: {{radius:0.15, colorscheme:'whiteCarbon'}} }});
            else if (style === 'stick') v.setStyle(proteinSel, {{stick: {{radius: 0.2, ...pStyle}}}});
            else if (style === 'surface') {{
                v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.8, ...pStyle}}, proteinSel, proteinSel);
                v.setStyle(proteinSel, {{cartoon: {{radius:0.1, color:'black', opacity:0.3}}}});
            }}
            
            // 2. Hotspots
            if (colorScheme === 'energy' && hotspots.length > 0) {{
                hotspots.forEach(function(h) {{
                   let color = 'white';
                   if (h.energy < -5.0) color = '#8B0000'; 
                   else if (h.energy < -2.0) color = 'red';
                   else if (h.energy < -1.0) color = 'orange';
                   else if (h.energy > 1.0) color = 'blue';
                   else return;
                   
                   let sel = {{resi: h.resi}};
                   
                   if(style !== 'stick') v.setStyle(sel, {{cartoon: {{color: color}}, stick: {{radius:0.2, color: color}} }});
                   else v.setStyle(sel, {{stick: {{radius:0.3, color: color}} }});
                   
                   if (showLabels && h.energy < -1.5) {{
                        v.addLabel(h.resn + h.resi + "\\n" + h.energy.toFixed(1), {{
                            position: sel, fontSize: 10, backgroundColor: color, fontColor: 'white', showBackground:true
                        }});
                   }}
                   
                   // 3. Energetic Links (Dashed Lines from Ligand to Residue)
                   if (showLinks && showLigand) {{
                       // We need centroid of residue and centroid of ligand
                       // Since we can't easily get centroids in this loop without async atoms,
                       // we'll use a simpler selection-based approach or addCylinder between atoms?
                       // Actually, we can just use "v.addCylinder" between selection centers
                       v.addCylinder({{
                           start: {{resi: h.resi, atom: 'CA'}}, // Alpha Carbon
                           end: {{ ...ligandSel }}, // Ligand Center (approx)
                           radius: 0.05,
                           color: color,
                           dashed: true,
                           fromCap: 1, toCap: 1
                       }});
                   }}
                }});
            }}
            
            // 4. Structural H-Bonds
            if (showHbonds && showLigand) {{
                if (typeof v.addHydrogenBonds === 'function') {{
                    try {{
                        v.addHydrogenBonds({{
                            sel0: ligandSel, sel1: proteinSel,
                            dist: 3.5, color: '#ffff00', dashed: false, linewidth: 2.0, alpha: 0.6
                        }});
                    }} catch (e) {{
                        console.warn("H-Bonds failed:", e);
                    }}
                }} else {{
                    console.warn("v.addHydrogenBonds is not supported in this 3Dmol version.");
                }}
            }}
            
            // 5. Ligand
            if (showLigand) {{
                v.addStyle(ligandSel, {{stick: {{colorscheme: 'greenCarbon', radius: 0.3}}}});
                if(showLabels) v.addLabel(ligVal || "LIG", {{position: ligandSel, backgroundColor:'black', fontColor:'green'}});
            }}
            
            v.render();
        }};
        
        // Tab Handling Logic (Vanilla JS)
        document.querySelectorAll('.nav-link').forEach(link => {{
            link.addEventListener('click', function(e) {{
                e.preventDefault();
                // Remove active class
                document.querySelectorAll('.nav-link').forEach(l => l.classList.remove('active'));
                document.querySelectorAll('.tab-pane').forEach(p => {{
                    p.classList.remove('show', 'active');
                    p.style.display = 'none';
                }});
                
                // Add active
                this.classList.add('active');
                let target = this.getAttribute('href');
                let pane = document.querySelector(target);
                pane.classList.add('show', 'active');
                pane.style.display = 'block';
            }});
        }});
        // Init state
        document.querySelector('#table').style.display = 'none'; // Hide table initially
        
        // Cleanup old
        window.takeSnapshot = function() {{
            window.viewer.render();
            let canvas = document.querySelector('#structure-viewer canvas');
            let img = canvas.toDataURL("image/png");
            let link = document.createElement('a');
            link.download = `mmgbsa_report_snap.png`;
            link.href = img;
            link.click();
        }};
        
        window.viewer.zoomTo();
        window.updateViz();
        """
        
        return script, controls_html

    def _generate_plotly_script(self, frame_data, analysis_results):
        """Generate Plotly.js code for the Heatmap."""
        if not frame_data:
            return ""
            
        try:
            df = pd.DataFrame(frame_data)
            
            # Identical logic to decomposition.py for picking top residues
            energy_cols = [col for col in df.columns if col.endswith('_total')]
            mean_energies = {col: df[col].mean() for col in energy_cols}
            
            # Sort and take top 20 (or matched to config if passed)
            top_n = 20
            sorted_residues = sorted(mean_energies.items(), key=lambda x: x[1])
            top_residues = sorted_residues[:top_n]
            
            # Add Ligand Total
            total_energy_series = df[energy_cols].sum(axis=1)
            total_mean = total_energy_series.mean()
            
            # Prepare Data Z-Matrix (Rows: Residues, Cols: Frames)
            z_data = []
            y_labels = []
            x_labels = df['frame_number'].tolist()
            
            for col, mean_val in top_residues:
                clean_name = col.replace('_total', '').replace('_A', '')
                y_labels.append(f"{clean_name} ({mean_val:.1f})")
                z_data.append(df[col].tolist())
            
            # Add Ligand
            z_data.append(total_energy_series.tolist())
            y_labels.append(f"LIGAND Total ({total_mean:.1f})")
            
            # Sanitize Z-Data (Replace NaNs with None for valid JSON)
            # Recursively or simply via list comprehension if needed, but since it's a list of lists of floats:
            z_data_sanitized = [[None if (val is None or np.isnan(val)) else val for val in row] for row in z_data]
            
            # Serialize for JS
            z_json = json.dumps(z_data_sanitized)
            x_json = json.dumps(x_labels)
            y_json = json.dumps(y_labels)
            
            script = f"""
            var data = [
                {{
                    z: {z_json},
                    x: {x_json},
                    y: {y_json},
                    type: 'heatmap',
                    colorscale: 'RdBu',
                    reversescale: true,
                    zmid: 0,
                    hoverongaps: false,
                    colorbar: {{
                        title: 'Energy (kcal/mol)',
                        titleside: 'right'
                    }}
                }}
            ];

            var layout = {{
                title: 'Per-Residue Interaction Energy Timeline',
                xaxis: {{ title: 'Frame Number' }},
                yaxis: {{ title: 'Residue', automargin: true }},
                margin: {{ l: 150, r: 50, b: 50, t: 50 }} 
            }};

            Plotly.newPlot('heatmap-plot', data, layout, {{responsive: true}});
            """
            return script
            
        except Exception as e:
            return f"// Failed to generate Plotly script: {e}"
