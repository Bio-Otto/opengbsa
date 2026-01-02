
import numpy as np
import mdtraj as md
from scipy.constants import k as k_B, h, N_A, c, calorie
from openmm import unit

class QuasiHarmonicAnalysis:
    """
    Perform Quasi-Harmonic Analysis (QHA) to calculate vibrational entropy
    from a molecular dynamics trajectory.
    """

    def __init__(self, trajectory, system, temperature=300.0, verbose=True):
        """
        Initialize QHA.

        Args:
            trajectory (mdtraj.Trajectory): Aligned trajectory of the solute (e.g., ligand).
            system (openmm.System): OpenMM system object (to get masses).
            temperature (float): Temperature in Kelvin.
        """
        self.trajectory = trajectory
        self.system = system
        self.temperature = temperature
        self.verbose = verbose
        
        # Physical Constants
        self.kB_kcal = k_B * N_A / (calorie * 1000)  # kcal/(mol*K) ≈ 0.001987
        self.h_kcal_s = h * N_A / (calorie * 1000)   # kcal*s/mol
        self.c_cm_s = c * 100                # cm/s

    def calculate_entropy(self):
        """
        Calculate quasi-harmonic entropy using mass-weighted covariance matrix.
        """
        if self.verbose:
            print("  • [QHA] Aligning trajectory...")
        
        # 1. Align Trajectory (Superposition)
        # Use the first frame as reference
        self.trajectory.superpose(self.trajectory, frame=0)
        
        if self.verbose:
            print(f"  • [QHA] Calculating covariance matrix for {self.trajectory.n_atoms} atoms over {self.trajectory.n_frames} frames...")

        # 2. Prepare Coordinate Array (Frames, 3N)
        # MDTraj: (n_frames, n_atoms, 3) in nm
        xyz = self.trajectory.xyz
        n_frames, n_atoms, _ = xyz.shape
        
        # 3. Get Masses (Daltons)
        masses = []
        for i in range(n_atoms):
            masses.append(self.system.getParticleMass(i).value_in_unit(unit.dalton))
        masses = np.array(masses)
        
        # 4. Mass-Weighting coordinates
        # q_i = x_i * sqrt(m_i)
        # We need to flatten to (n_frames, 3*n_atoms)
        
        # Repeat masses for x, y, z
        mass_weights = np.repeat(np.sqrt(masses), 3) # Shape (3N,)
        
        # Flatten coordinates: (Frame, Atom, XYZ) -> (Frame, 3N)
        coords_flat = xyz.reshape(n_frames, -1)
        
        # Subtract mean position (optional if np.cov handles it, but good for mass-weighting center)
        mean_coords = np.mean(coords_flat, axis=0)
        diff_coords = coords_flat - mean_coords
        
        # Apply mass weights
        # mw_coords[t, i] = diff_coords[t, i] * mass_weights[i]
        mw_coords = diff_coords * mass_weights
        
        # 5. Calculate Covariance Matrix (3N x 3N)
        # Units: (nm * sqrt(Dalton))^2
        sigma = np.cov(mw_coords, rowvar=False)
        
        # 6. Diagonalize
        if self.verbose:
            print("  • [QHA] Diagonalizing covariance matrix...")
        
        eigenvalues, _ = np.linalg.eigh(sigma)
        
        # Sort eigenvalues (largest variance = lowest frequency)
        # We need to check units carefully.
        # Variance sigma_i = k_B * T / omega_i^2
        # So omega_i = sqrt(k_B * T / sigma_i)
        
        # Unit conversion: 
        # Sigma units: nm^2 * Dalton
        # k_B * T units: J/molecule? No, we work in consistent units.
        
        # Let's perform calculation in SI units first to be safe.
        # k_B (SI) = 1.380649e-23 J/K
        # T (K) = 300
        # Sigma (nm^2 * Da) -> Convert to m^2 * kg
        # 1 nm = 1e-9 m
        # 1 Da = 1.66053907e-27 kg
        
        factor_dist = 1e-9
        factor_mass = 1.66053907e-27
        
        sigma_SI = eigenvalues * (factor_dist**2) * factor_mass
        
        k_B_SI = k_B # J/K
        kbT = k_B_SI * self.temperature
        
        # Calculate Angular Frequencies (omega)
        # omega = sqrt(kbT / lambda)
        # Handle small/negative eigenvalues (rot/trans modes + numerical noise)
        
        # Ideally, we should remove 6 smallest FREQUENCIES (which correspond to LARGEST variances/eigenvalues?)
        # Wait. 
        # Large variance -> Soft mode -> Low frequency.
        # Small variance -> Stiff mode -> High frequency.
        # Zero variance -> Infinite frequency? No.
        
        # Rotational/Translational modes have HUGE variance (infinite in unbound, but we aligned).
        # Alignment removes rigid body motion, so the 6 rot/trans modes should have near-zero eigenvalues (variance).
        # Wait, if we remove motion, the variance is ZERO. 
        # So the eigenvalues for rot/trans should be approx ZERO.
        # But wait, omega^2 = kT/lambda. If lambda -> 0, omega -> Infinity (stiff?).
        
        # Let's re-read QHA theory.
        # "The eigenvectors with the largest eigenvalues correspond to the large-amplitude motions (low frequency)."
        # "The eigenvectors with the smallest eigenvalues correspond to high-frequency fluctuations."
        # Rotational/Translational modes, if not removed, would have infinite variance.
        # Since we aligned, they should have near 0 variance? 
        # Actually, superposition minimizes RMSD, effectively removing the "variance" associated with rotation/translation.
        # So we expect 6 eigenvalues to be very close to zero.
        # omega = sqrt(kT/lambda).
        # If lambda is tiny, omega is huge? That doesn't sound like "zero frequency mode".
        
        # Standard QHA formulation (Karplus & Kushick, 1981): 
        # S = k/2 * ln(det(1 + ...)) -> approximation using eigenvalues.
        # S_vib = sum( S_ho(omega_i) )
        # where omega_i = sqrt(k_B T / lambda_i).
        
        # If lambda_i (variance) is small -> omega is large (stiff bond).
        # If lambda_i is large -> omega is small (soft mode).
        
        # What about the 6 zero modes?
        # Superposition creates 6 degrees of freedom with 0 variance (eigenvalues ~ 0).
        # sqrt(kT / 0) -> Infinity.
        # Infinite frequency -> Entropy = 0.
        # Correct. We should exclude them or let them be 0 entropy.
        
        # So we process all non-zero eigenvalues.
        
        entropies = []
        frequencies_cm = []
        
        # Cutoff logic (optional but recommended to avoid noise)
        # Similar to NMA, we can ignore very small eigenvalues (which mean infinite freq -> 0 entropy).
        
        if self.verbose:
             # print(f"  • [QHA] Computed TS = {TS:.4f} kcal/mol")
             pass
        
        for val in sigma_SI:
            if val < 1e-65: # Lowered threshold to account for SI scaling
                continue
            
            omega = np.sqrt(kbT / val) # rad/s
            
            # Convert to cm^-1 for reporting
            # v (Hz) = omega / 2pi
            # wavenumber = v / c
            freq_cm = (omega / (2 * np.pi * c * 100))
            frequencies_cm.append(freq_cm)
            
            # Calculate Entropy for this mode using Quantum Harmonic Oscillator formula
            # x = h_bar * omega / (k_B * T)
            # h_bar = h / 2pi
            
            hbar = h / (2 * np.pi)
            x = (hbar * omega) / (k_B_SI * self.temperature)
            
            # x = E / kT
            # S = R * [ x/(e^x - 1) - ln(1 - e^-x) ] (per mole if R) or k if per particle
            
            if x > 50: # High frequency limit -> S ~ 0
                term = 0.0
            else:
                try:
                    term = (x / (np.exp(x) - 1.0)) - np.log(1.0 - np.exp(-x))
                except:
                    term = 0.0
            
            entropies.append(term)
            
        total_entropy_cal_mol_k = np.sum(entropies) * self.kB_kcal * 1000 # Convert k_B to kcal/mol/K ?? 
        # self.kB_kcal is kcal/(mol*K). Result of sum is dimensionless factor (S/k).
        # So multiply by (R in kcal/mol/K).
        # self.kB_kcal is ~0.001987.
        
        S_total = np.sum(entropies) * self.kB_kcal
        
        # T*S in kcal/mol
        TS = self.temperature * S_total
        
        if self.verbose:
             print(f"  • [QHA] Computed TS = {TS:.4f} kcal/mol")
             
        return TS

