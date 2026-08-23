"""
Binding free energy -> Kd/IC50/pIC50 conversion helpers.

See `calculate_ic50`'s docstring for an important caveat: despite the
name, it reports a thermodynamic Kd (via Kd = exp(dG/RT)), not a true
Cheng-Prusoff-corrected IC50, and is not currently wired into the main
`mmgbsa.runner`/`mmgbsa.cli` pipeline.
"""
import numpy as np
from openmm import unit

def calculate_ic50(delta_g, delta_g_std, temperature=300.0*unit.kelvin):
    """
    Calculate a theoretical Kd (reported as 'ic50') from Binding Free Energy,
    with a confidence interval from `delta_g_std`.

    NOTE: despite the name, this computes Kd = exp(delta_g / RT), the
    thermodynamic dissociation constant -- NOT a true IC50. For a competitive
    inhibitor, IC50 relates to Kd via the Cheng-Prusoff equation,
    IC50 = Kd * (1 + [S]/Km), which depends on substrate concentration and Km
    and is NOT applied here. Treat this function's output as a Kd estimate,
    not a literal IC50 prediction, unless [S]/Km happens to be negligible.

    Also note: as of this writing, neither this function nor the module that
    calls it (`mmgbsa.plotting`) is imported/reachable from the main
    pipeline (`mmgbsa.runner`/`mmgbsa.cli`) -- this is effectively unused
    code. Verify it is wired in before relying on its output.

    Parameters
    ----------
    delta_g : float or unit.Quantity
        Mean binding free energy (kcal/mol)
    delta_g_std : float or unit.Quantity
        Standard deviation or Standard Error of the Mean (kcal/mol)
    temperature : unit.Quantity
        Temperature (default: 300K)

    Returns
    -------
    dict
        Dictionary containing:
        - 'ic50': Theoretical Kd, labeled IC50 (micromolar)
        - 'ic50_low': Lower bound of 95% CI
        - 'ic50_high': Upper bound of 95% CI
        - 'unit': 'micromolar'
    """
    
    # Constants
    R = 1.98720425864083e-3 # kcal/(mol*K)
    
    # Ensure units are handled
    if unit.is_quantity(delta_g):
        delta_g = delta_g.value_in_unit(unit.kilocalories_per_mole)
    if unit.is_quantity(delta_g_std):
        delta_g_std = delta_g_std.value_in_unit(unit.kilocalories_per_mole)
    if unit.is_quantity(temperature):
        T = temperature.value_in_unit(unit.kelvin)
    else:
        T = temperature
        
    RT = R * T
    
    # Z-score for 95% CI
    z = 1.96
    
    # Calculate bounds in Energy space first
    # Note: Higher Energy = Weaker Binding = Higher IC50
    # So Low Energy (Strong) -> Low IC50
    # High Energy (Weak) -> High IC50
    
    dg_low = delta_g - (z * delta_g_std)
    dg_high = delta_g + (z * delta_g_std)
    
    # Calculate IC50 = exp(dG / RT) * 1M -> convert to uM (1e6)
    # The relation is Kd = exp(dG/RT). For inhibitors, IC50 ~ Kd (Cheng-Prusoff for competitive)
    
    def to_ic50_uM(dg_val):
        """Convert a binding free energy (kcal/mol) to a Kd estimate in uM."""
        # Kd in Molar
        kd_molar = np.exp(dg_val / RT)
        # Convert to micromolar
        return kd_molar * 1e6
        
    ic50_mean = to_ic50_uM(delta_g)
    ic50_lower_bound = to_ic50_uM(dg_low) # Lower dG -> Lower IC50
    ic50_upper_bound = to_ic50_uM(dg_high) # Higher dG -> Higher IC50
    
    
    # Calculate pIC50 = -log10(IC50_Molar)
    # IC50_Molar = exp(dG / RT)
    # pIC50 = - (dG / RT) / ln(10)
    #       = - dG / (RT * 2.303)
    
    pic50_mean = -delta_g / (RT * np.log(10))
    pic50_low = -dg_high / (RT * np.log(10)) # High dG (weak) -> Low pIC50
    pic50_high = -dg_low / (RT * np.log(10)) # Low dG (strong) -> High pIC50
    
    return {
        'ic50': ic50_mean,
        'ic50_low': ic50_lower_bound,
        'ic50_high': ic50_upper_bound,
        'unit': 'micromolar',
        'pic50': pic50_mean,
        'pic50_low': pic50_low,
        'pic50_high': pic50_high
    }
