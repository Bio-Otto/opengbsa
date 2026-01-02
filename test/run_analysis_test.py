#!/usr/bin/env python3
import sys
import yaml
import os
from mmgbsa.runner import MMGBSARunner
from mmgbsa.logger import ToolLogger

log = ToolLogger()

def main():
    if len(sys.argv) < 2:
        log.error("Usage: python run_analysis_test.py <config.yaml>")
        sys.exit(1)
        
    config_file = sys.argv[1]
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
        
    runner = MMGBSARunner(config)
    results = runner.run_analysis()
    
    if results:
        log.section("FINAL RESULTS SUMMARY")
        if 'mean_binding_energy' in results:
            log.result("Binding Energy", f"{results['mean_binding_energy']:.2f} ± {results['std_error']:.2f}", "kcal/mol")
        
        if 'entropy_term' in results:
            log.result("Entropy Term", f"{results['entropy_term']:.2f}", "kcal/mol")
            
        if 'delta_g_total' in results:
            log.result("Total Free Energy", f"{results['delta_g_total']:.2f}", "kcal/mol")

if __name__ == "__main__":
    main()
