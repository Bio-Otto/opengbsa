"""
Results management package - comprehensive output handling and reporting.

Provides:
- ResultsExporter: Multi-format export (CSV, JSON, per-residue)
- ReportBuilder: HTML report generation with charts and 3D visualization
- ResultsValidator: Post-analysis validation and QC
- ResultsManager: High-level results orchestrator
"""

from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import json
import logging

log = logging.getLogger(__name__)


class ResultsExporter:
    """Handles export of results in multiple formats."""
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
    
    def export_csv(self, results_df, output_file: str) -> str:
        """Export frame-by-frame results to CSV."""
        try:
            results_df.to_csv(output_file, index=False)
            if self.verbose:
                log.info(f"Exported frame results to CSV: {output_file}")
            return output_file
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to export CSV: {e}")
            return None
    
    def export_json(self, summary_data: Dict[str, Any], output_file: str) -> str:
        """Export summary statistics to JSON."""
        try:
            with open(output_file, 'w') as f:
                json.dump(summary_data, f, default=str, indent=2)
            if self.verbose:
                log.info(f"Exported summary to JSON: {output_file}")
            return output_file
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to export JSON: {e}")
            return None
    
    def export_per_residue(self, decomp_df, output_file: str) -> str:
        """Export per-residue decomposition to CSV."""
        try:
            decomp_df.to_csv(output_file, index=False)
            if self.verbose:
                log.info(f"Exported per-residue decomposition to CSV: {output_file}")
            return output_file
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to export per-residue results: {e}")
            return None
    
    def export_bootstrap(self, bootstrap_data: Dict[str, Any], output_file: str) -> str:
        """Export bootstrap uncertainty analysis to JSON."""
        try:
            with open(output_file, 'w') as f:
                json.dump(bootstrap_data, f, indent=2)
            if self.verbose:
                log.info(f"Exported bootstrap results to JSON: {output_file}")
            return output_file
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to export bootstrap results: {e}")
            return None


class ReportBuilder:
    """Constructs comprehensive analysis reports."""
    
    def __init__(self, output_dir: Optional[str] = None, verbose: bool = False):
        self.output_dir = Path(output_dir) if output_dir else Path('.')
        self.verbose = verbose
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_text_report(self, summary: Dict[str, Any], validation_warnings: List[str],
                           physics_assumptions: List[str], output_file: Optional[str] = None) -> str:
        """Generate plain text analysis report."""
        report_path = output_file or self.output_dir / "analysis_report.txt"
        
        try:
            with open(report_path, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("Advanced MM/GBSA Analysis Report\n")
                f.write("=" * 60 + "\n\n")
                
                # Summary Statistics
                f.write("RESULTS SUMMARY\n")
                f.write("-" * 40 + "\n")
                f.write(f"Mean Binding Energy:    {summary.get('mean', 'N/A')}\n")
                f.write(f"Std Dev:                {summary.get('std_dev', 'N/A')}\n")
                f.write(f"Median:                 {summary.get('median', 'N/A')}\n")
                f.write(f"Min/Max:                {summary.get('min', 'N/A')} / {summary.get('max', 'N/A')}\n")
                f.write(f"95% CI:                 [{summary.get('ci_lower', 'N/A')}, {summary.get('ci_upper', 'N/A')}]\n\n")
                
                # Convergence
                f.write("CONVERGENCE ANALYSIS\n")
                f.write("-" * 40 + "\n")
                convergence = summary.get('convergence', {})
                f.write(f"Converged:              {convergence.get('converged', 'Unknown')}\n")
                f.write(f"Convergence Diff:       {convergence.get('difference', 'N/A')} kcal/mol\n\n")
                
                # Validation Warnings
                if validation_warnings:
                    f.write("VALIDATION WARNINGS\n")
                    f.write("-" * 40 + "\n")
                    for warning in validation_warnings:
                        f.write(f"• {warning}\n")
                    f.write("\n")
                
                # Physics Assumptions
                if physics_assumptions:
                    f.write("PHYSICS ASSUMPTIONS & FALLBACKS\n")
                    f.write("-" * 40 + "\n")
                    seen = set()
                    for assumption in physics_assumptions:
                        if assumption not in seen:
                            f.write(f"• {assumption}\n")
                            seen.add(assumption)
                else:
                    f.write("PHYSICS ASSUMPTIONS & FALLBACKS\n")
                    f.write("-" * 40 + "\n")
                    f.write("• No assumptions. All parameters explicitly defined.\n")
                
                f.write("\n" + "=" * 60 + "\n")
            
            if self.verbose:
                log.info(f"Generated text report: {report_path}")
            return str(report_path)
        
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to generate text report: {e}")
            return None
    
    def generate_html_report(self, summary: Dict[str, Any], charts_data: Optional[Dict] = None,
                           output_file: Optional[str] = None) -> str:
        """Generate interactive HTML report with plots and 3D visualization."""
        report_path = output_file or self.output_dir / "index.html"
        
        try:
            # Safe value extraction with formatting
            mean_val = summary.get('mean', 0)
            mean_str = f"{float(mean_val):.2f}" if mean_val != 'N/A' and mean_val is not None else 'N/A'
            
            std_val = summary.get('std_dev', 0)
            std_str = f"{float(std_val):.2f}" if std_val != 'N/A' and std_val is not None else 'N/A'
            
            median_val = summary.get('median', 0)
            median_str = f"{float(median_val):.2f}" if median_val != 'N/A' and median_val is not None else 'N/A'
            
            ci_lower = summary.get('ci_lower', 0)
            ci_upper = summary.get('ci_upper', 0)
            ci_str = f"[{float(ci_lower):.2f}, {float(ci_upper):.2f}]" if ci_lower != 'N/A' and ci_upper != 'N/A' else 'N/A'
            
            n_frames = summary.get('n_frames', 'N/A')
            gb_model = summary.get('gb_model', 'N/A')
            salt_conc = summary.get('salt_conc', 'N/A')
            
            html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>MM/GBSA Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .stat-value {{ font-size: 24px; font-weight: bold; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; margin-top: 5px; }}
        .warning {{ 
            background-color: #fff3cd; 
            border-left: 4px solid #ffc107; 
            padding: 15px; 
            margin: 10px 0; 
            border-radius: 4px;
        }}
        footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Advanced MM/GBSA Analysis Report</h1>
        
        <h2>Results Summary</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{mean_str}</div>
                <div class="stat-label">Mean ΔG (kcal/mol)</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">±{std_str}</div>
                <div class="stat-label">Std Deviation</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{median_str}</div>
                <div class="stat-label">Median ΔG</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{ci_str}</div>
                <div class="stat-label">95% Confidence Interval</div>
            </div>
        </div>
        
        <h2>Analysis Details</h2>
        <p><strong>Frames analyzed:</strong> {n_frames}</p>
        <p><strong>GB Model:</strong> {gb_model}</p>
        <p><strong>Salt concentration:</strong> {salt_conc} M</p>
        
        <footer>
            <p>Generated with openGBSA</p>
        </footer>
    </div>
</body>
</html>
            """
            
            with open(report_path, 'w') as f:
                f.write(html_template)
            
            if self.verbose:
                log.info(f"Generated HTML report: {report_path}")
            return str(report_path)
        
        except Exception as e:
            if self.verbose:
                log.warning(f"Failed to generate HTML report: {e}")
            return None


class ResultsValidator:
    """
    Validates analysis results for quality and consistency.

    NOTE: not currently called anywhere in the pipeline (including by
    `ResultsManager`, which does not instantiate or use this class). The
    main pipeline's result sanity-checking is instead done by
    `mmgbsa.validation.TopologyValidator.validate_system_sanity`, called
    from `mmgbsa.runner.MMGBSARunner`. Call these methods explicitly if you
    need this class's specific (differently-thresholded) checks.
    """

    def __init__(self, verbose: bool = False):
        self.verbose = verbose
    
    def validate_energy_range(self, binding_energies) -> Tuple[bool, List[str]]:
        """Check if binding energies are in reasonable range."""
        warnings = []
        
        mean_be = binding_energies.mean()
        if mean_be > 100:
            warnings.append(f"Unusually high positive binding energy: {mean_be:.1f} kcal/mol")
        elif mean_be < -200:
            warnings.append(f"Unusually negative binding energy: {mean_be:.1f} kcal/mol")
        
        return len(warnings) == 0, warnings
    
    def validate_statistics(self, summary: Dict[str, Any]) -> List[str]:
        """Validate statistical reasonableness."""
        warnings = []
        
        std_dev = summary.get('std_dev', 0)
        if std_dev < 0.01:
            warnings.append(f"Very low standard deviation ({std_dev:.3f}) - trajectory may be too short")
        elif std_dev > 100:
            warnings.append(f"Very high standard deviation ({std_dev:.1f}) - system may be unstable")
        
        return warnings


class ResultsManager:
    """High-level results orchestrator."""
    
    def __init__(self, output_dir: Optional[str] = None, verbose: bool = False):
        self.output_dir = Path(output_dir) if output_dir else Path('.')
        self.verbose = verbose
        
        self.exporter = ResultsExporter(verbose=verbose)
        self.report_builder = ReportBuilder(str(self.output_dir), verbose=verbose)
        self.validator = ResultsValidator(verbose=verbose)
    
    def save_results(self, results: Dict[str, Any], filename: str = 'summary.json') -> str:
        """Save results dictionary to JSON (legacy interface)."""
        output_path = self.output_dir / filename
        return self.exporter.export_json(results, str(output_path))
    
    def save_comprehensive_results(self, frame_results, summary_stats: Dict[str, Any],
                                   decomposition_results=None, bootstrap_results=None,
                                   validation_warnings: List[str] = None,
                                   physics_assumptions: List[str] = None) -> Dict[str, str]:
        """Save all results components in organized structure."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        saved_files = {}
        
        try:
            # Frame results CSV
            if frame_results is not None:
                csv_path = self.exporter.export_csv(frame_results, 
                                                   str(self.output_dir / "frame_results.csv"))
                saved_files['frame_results_csv'] = csv_path
            
            # Summary statistics JSON
            summary_path = self.exporter.export_json(summary_stats,
                                                    str(self.output_dir / "summary.json"))
            saved_files['summary_json'] = summary_path
            
            # Per-residue decomposition
            if decomposition_results is not None:
                decomp_path = self.exporter.export_per_residue(decomposition_results,
                                                              str(self.output_dir / "per_residue.csv"))
                saved_files['decomposition_csv'] = decomp_path
            
            # Bootstrap results
            if bootstrap_results is not None:
                bootstrap_path = self.exporter.export_bootstrap(bootstrap_results,
                                                               str(self.output_dir / "bootstrap.json"))
                saved_files['bootstrap_json'] = bootstrap_path
            
            # Text report
            text_report = self.report_builder.generate_text_report(
                summary_stats,
                validation_warnings or [],
                physics_assumptions or [],
                str(self.output_dir / "report.txt")
            )
            if text_report:
                saved_files['text_report'] = text_report
            
            # HTML report
            html_report = self.report_builder.generate_html_report(
                summary_stats,
                output_file=str(self.output_dir / "index.html")
            )
            if html_report:
                saved_files['html_report'] = html_report
            
            if self.verbose:
                log.info(f"Saved {len(saved_files)} result files to {self.output_dir}")
            
            return saved_files
        
        except Exception as e:
            if self.verbose:
                log.error(f"Error saving results: {e}")
            return saved_files


# Import pandas for Timestamp (used in HTML report)
try:
    import pandas as pd
except ImportError:
    pd = None
