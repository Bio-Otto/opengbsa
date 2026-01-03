
import pandas as pd
import numpy as np
from mmgbsa.core import GBSACalculator

def test_entropy():
    print("Loading results...")
    try:
        df = pd.read_csv("fixed_enhanced_mmgbsa_results_gbn2.csv")
        binding_energies = df['binding_energy'].values
        print(f"Binding energies: {binding_energies}")
        print(f"Shape: {binding_energies.shape}")
        
        calculator = GBSACalculator(temperature=300)
        
        print("Calculating entropy...")
        entropy = calculator.calculate_interaction_entropy(binding_energies, 300.0)
        print(f"Entropy: {entropy}")
        
    except Exception as e:
        print(f"Caught expected error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_entropy()
