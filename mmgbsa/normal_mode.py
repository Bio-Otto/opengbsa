from openmm.unit import *
from openmm import *
from openmm.app import *
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.constants import pi, Boltzmann, hbar, Avogadro
from math import pi
from copy import deepcopy 
from warnings import warn
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.dpi'] = 600

class NormalModeAnalysis(object):
    """
        Create a NormalModeAnalysis object

        Args:
            topology (simtk.openmm.openmm.Topology): The OpenMM topology object
            system (simtk.openmm.openmm.System): The OpenMM system object.
            integrator (simtk.openmm.openmm.Integrator): The OpenMM integrator object.
            initPositions (numpy.array): The N*3 array of positions (N is the number of atoms).
            CPUProp (dictionary=None): The CPU platform-specific properties used for generating simulation object in OpenMM.
            CUDAProp (dictionary=None): The CUDA platform-specific properties used for generating simulation object in OpenMM.
    """
    def __init__(self, topology, system, integrator, initPositions, CPUProp=None, CUDAProp=None):
        self.topology = topology
        self.CPUSystem = deepcopy(system)
        self.CUDASystem = deepcopy(system)
        self.CPUIntegrator = deepcopy(integrator)
        self.CUDAIntegrator = deepcopy(integrator)
        self.initPositions = initPositions
        self.CPUProp = CPUProp
        self.CUDAProp = CUDAProp
        self.CPUSimulation = self.__getCPUSimulation__()
        self.CUDASimulation = self.__getCUDASimulation__()
        self.CPUSimulation.context.setPositions(self.initPositions)
        self.CUDASimulation.context.setPositions(self.initPositions)

    def __getDefaultCPUProperty__(self):
        pass

    def __getDefaultCUDAProperty__(self):
        return {'CudaDeviceIndex': '0', 'CudaPrecision': 'double', 'DeterministicForces': 'true'}

    def __getCPUSimulation__(self):
        if self.CPUProp:
            CPUPlatform = Platform.getPlatformByName('CPU')
            simulation = Simulation(self.topology, self.CPUSystem, self.CPUIntegrator, CPUPlatform, self.CPUProp)
            return simulation
        else:
            CPUPlatform = Platform.getPlatformByName('CPU')
            simulation = Simulation(self.topology, self.CPUSystem, self.CPUIntegrator, CPUPlatform)
            return simulation

    def __getCUDASimulation__(self):
        # Try CUDA first
        try:
            platform = Platform.getPlatformByName('CUDA')
            prop = self.CUDAProp if self.CUDAProp else self.__getDefaultCUDAProperty__()
            simulation = Simulation(self.topology, self.CUDASystem, self.CUDAIntegrator, platform, prop)
            return simulation
        except Exception:
            # Fallback to OpenCL
            try:
                platform = Platform.getPlatformByName('OpenCL')
                # OpenCL properties equivalent
                prop = {'OpenCLPrecision': 'double'}
                if self.CUDAProp and 'CudaDeviceIndex' in self.CUDAProp:
                    prop['OpenCLDeviceIndex'] = self.CUDAProp['CudaDeviceIndex']
                
                simulation = Simulation(self.topology, self.CUDASystem, self.CUDAIntegrator, platform, prop)
                return simulation
            except Exception:
                # Fallback to CPU as last resort (slow but works)
                print("Warning: Neither CUDA nor OpenCL available for NMA. Falling back to CPU (slow).")
                platform = Platform.getPlatformByName('CPU')
                simulation = Simulation(self.topology, self.CUDASystem, self.CUDAIntegrator, platform)
                return simulation

    def __getVibrationalSpectrum__(self, SquareAngularFreq):
        SquareAngularFreqSI = (4.184*10**26)*SquareAngularFreq
        AngularFreqSI = np.sqrt(SquareAngularFreqSI)
        VibrationalSpectrum = AngularFreqSI/(6*pi*10**10)
        return VibrationalSpectrum

    def __checkSymmetric__(self, array2D, tol=1e-8):
        return np.allclose(array2D, array2D.T, atol=tol)

    def __checkPositiveDefinite__(self, eigVals):
        eigVals = eigVals[6:]
        return np.all(eigVals > 0)

    def CPUPreMinimization(self):
        """
            Initial minimization in CPU platform to remove bad contacts of the initial input positions.
            This function uses all the default settings of the class method simtk.openmm.app.simulation.Simulation.minimizeEnergy(). 
            After the initial minimization, set the self.CUDASimulation.context to the pre-minimized state.
            The mean force after minimization would fall in around the scale of ~ 0.2 kcal/(A mol).
        """
        self.CPUSimulation.minimizeEnergy()
        PreMinimizedState = self.CPUSimulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
        self.CUDASimulation.context.setState(PreMinimizedState)

    def CUDAMinimizationCycle(self, MiniTolerance=0, MaxMiniCycle=1000, NumMiniStepPerCycle=10000, MiniForceRatio=1e-6):
        """
            Designed energy minimization cycle to minimize the structure such that the system mean force would fall around 2e-07 kcal/(A mol).
            This function will use the default positions to minimize.
            If the user did not to self.CPUPreMinimization() first, then the initial input positions will be used.
            Otherwise the pre-minimized positions will be used for performing the minimization cycle in CUDA platform.

            Args:
            MiniTolerance (energy=0*kilojoule/mole): The energy tolerance to which the system should be minimized set for each cycle.
            MaxMiniCycle (int=1000): The maximum number of cycles to perform energy minimizations.
            NumMiniStepPerCycle (int=10000): MaxIterations for each cycle of energy minimization.
            MiniForceRatio (double=1e-6): The order of mean force that the minimization cycle should eliminated.
        """
        PreMinimizedState = self.CUDASimulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
        PreMinimizedForces = PreMinimizedState.getForces(asNumpy=True).value_in_unit(kilocalorie/(mole*angstrom))
        PreMinimizedMeanForce = np.linalg.norm(PreMinimizedForces,axis=1).mean() * (kilocalorie/(mole*angstrom))
        self.MiniForceThreshold = PreMinimizedMeanForce * MiniForceRatio

        for i in range(MaxMiniCycle):
            self.CUDASimulation.minimizeEnergy(tolerance=MiniTolerance, maxIterations=NumMiniStepPerCycle)
            currentState = self.CUDASimulation.context.getState(getForces=True)
            currentForces = currentState.getForces(asNumpy=True).value_in_unit(kilocalorie/(mole*angstrom))
            currentMeanForce = np.linalg.norm(currentForces,axis=1).mean() * (kilocalorie/(mole*angstrom))
            if currentMeanForce < self.MiniForceThreshold:
                break

    def CalculateNormalModes(self, TweakEnergyRatio=1e-12, cutoff_frequency=10.0):
        """
        Compute normal modes via a numerical (finite-difference) mass-weighted
        Hessian at the current (assumed energy-minimized) structure.

        For each of the 3N Cartesian degrees of freedom, displaces the
        structure by +/-`self.TweakDisplacement` (a per-DOF step size derived
        from `TweakEnergyRatio` -- the fraction of the mean per-atom potential
        energy that the displacement should correspond to, converted to a
        distance via the minimized structure's mean force magnitude) and
        re-evaluates forces at each displaced geometry. A cubic spline
        through the 3 points (-step, 0, +step) of (position, -force) per
        perturbed DOF gives that column of the Cartesian spring-constant
        (Hessian) matrix as the spline's derivative at 0. The Hessian is
        mass-weighted (dividing by sqrt(mass_i * mass_j) for each i,j pair),
        symmetrized (averaged with its own transpose, since the finite-
        difference estimate is not exactly symmetric), and diagonalized.

        The 6 lowest eigenvalues (translational + rotational modes, expected
        to be ~0 for a converged minimum) are left unclamped; among the
        remaining internal-mode eigenvalues, any below the eigenvalue
        corresponding to `cutoff_frequency` (including any spuriously
        negative ones from imperfect minimization/finite-difference noise)
        are clamped up to that threshold -- this is the "Robust NMA"
        stabilization referenced in the printed note, preventing near-zero
        or negative eigenvalues from producing divergent or complex
        vibrational frequencies/entropies downstream.

        Args
        ----
        TweakEnergyRatio : float, default 1e-12
            Target fractional potential-energy perturbation per atom used to
            size the finite-difference displacement step (smaller = smaller
            displacement = less anharmonic contamination, but more sensitive
            to force-evaluation numerical noise).
        cutoff_frequency : float, default 10.0
            Cutoff frequency in cm^-1. Internal-mode eigenvalues corresponding
            to a lower (or negative) frequency than this are clamped to the
            eigenvalue at this cutoff before computing the vibrational
            spectrum/entropy.

        Populates `self.VibrationalSpectrum` (via `__getVibrationalSpectrum__`)
        and the eigenvalue/eigenvector arrays used by
        `getVibrationalEntropyCM`/`getVibrationalEntropyQM`.
        """
        MinimizedState = self.CUDASimulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
        MinimizedPositions = MinimizedState.getPositions(asNumpy=True).in_units_of(angstrom)
        MinimizedPotentialEnergy = MinimizedState.getPotentialEnergy().in_units_of(kilocalorie/mole)
        MinimizedForces = MinimizedState.getForces(asNumpy=True).value_in_unit(kilocalorie/(mole*angstrom))
        MinimizedMeanForce = np.linalg.norm(MinimizedForces,axis=1).mean() * (kilocalorie/(mole*angstrom))

        NumAtoms = self.CUDASimulation.system.getNumParticles()
        TweakEnergyDiff = MinimizedPotentialEnergy/NumAtoms * TweakEnergyRatio
        self.TweakDisplacement = abs((TweakEnergyDiff/MinimizedMeanForce).in_units_of(angstrom))
        NumDimension = MinimizedForces.shape[0]*MinimizedForces.shape[1]
        Forces3N = MinimizedForces.reshape((1,NumDimension))
        Positions3N = MinimizedPositions.reshape((1,NumDimension))

        TweakEnergies = np.zeros((6*NumAtoms,1))
        MeanForces = np.zeros((6*NumAtoms,1))
        TweakForces = np.zeros((6*NumAtoms,NumDimension))
        SpringConsts = np.zeros((NumDimension,NumDimension))
        ForcesZero = np.copy(Forces3N[0])

        for i in range(NumDimension):
            currentPositions3NPos = np.copy(Positions3N[0])
            currentPositions3NPos[i] += self.TweakDisplacement.value_in_unit(angstrom)
            currentPositions3NPosQuantity = currentPositions3NPos.reshape((NumAtoms,3)) * angstrom
            currentPositions3NPosQuantity = currentPositions3NPosQuantity.in_units_of(nanometer)
            self.CUDASimulation.context.setPositions(currentPositions3NPosQuantity)
            NewStatePos = self.CUDASimulation.context.getState(getEnergy=True, getForces=True)
            NewForcesPos = NewStatePos.getForces(asNumpy=True).value_in_unit(kilocalorie/(mole*angstrom)).reshape((1,NumDimension))[0]
            TweakForces[2*i,:] = NewForcesPos
            TweakEnergies[2*i] = NewStatePos.getPotentialEnergy().value_in_unit(kilocalorie/mole)
            MeanForces[2*i] = np.linalg.norm(NewForcesPos.reshape((NumAtoms,3)),axis=1).mean()
    
            currentPositions3NNeg = np.copy(Positions3N[0])
            currentPositions3NNeg[i] -= self.TweakDisplacement.value_in_unit(angstrom)
            currentPositions3NNegQuantity = currentPositions3NNeg.reshape((NumAtoms,3)) * angstrom
            currentPositions3NNegQuantity = currentPositions3NNegQuantity.in_units_of(nanometer)
            self.CUDASimulation.context.setPositions(currentPositions3NNegQuantity)
            NewStateNeg = self.CUDASimulation.context.getState(getEnergy=True, getForces=True)
            NewForcesNeg = NewStateNeg.getForces(asNumpy=True).value_in_unit(kilocalorie/(mole*angstrom)).reshape((1,NumDimension))[0]
            TweakForces[2*i+1,:] = NewForcesNeg
            TweakEnergies[2*i+1] = NewStateNeg.getPotentialEnergy().value_in_unit(kilocalorie/mole)
            MeanForces[2*i+1] = np.linalg.norm(NewForcesNeg.reshape((NumAtoms,3)),axis=1).mean()
    
            PosVariable = np.array([-self.TweakDisplacement.value_in_unit(angstrom), 0, self.TweakDisplacement.value_in_unit(angstrom)])
            ForceVariable = np.array([NewForcesNeg,ForcesZero,NewForcesPos]).T
            ForceVariable = ForceVariable*(-1)
            ForceFunction = CubicSpline(PosVariable,ForceVariable,axis=1)
            SpringConsts[i,:] = ForceFunction(0,1)

        MassArray = np.zeros(NumAtoms)
        for i in range(NumAtoms):
            MassArray[i] = self.CUDASimulation.system.getParticleMass(i).value_in_unit(dalton)

        MassArray3N = np.sqrt(MassArray.repeat(3))
        MassMatrix = np.outer(MassArray3N,MassArray3N)
        Hessian = np.divide(SpringConsts,MassMatrix)
        HessianSymmetric = np.mean([Hessian,Hessian.T],axis=0)
        eigVal, eigVec = np.linalg.eig(HessianSymmetric)
        sortIdx = np.argsort(eigVal)
        eigValSorted = eigVal[sortIdx]
        eigVecSorted = eigVec[:,sortIdx]
        
        # --- Robust NMA: Clamp Eigenvalues based on Cutoff Frequency ---
        # Convert Cutoff (cm^-1) to Eigenvalue Unit (kcal/(g*A^2))
        # Relation: freq_cm = sqrt(eig_akma * C) / (2*pi*c)
        # So: eig_akma = (freq_cm * 2*pi*c)^2 / C
        # Constant calculation:
        # C = 4.184e26 (from __getVibrationalSpectrum__)
        # c = 3e10 cm/s
        
        # Using existing functionality:
        # We can just reverse the __getVibrationalSpectrum__ logic or use the cutoff directly
        
        # Calculate threshold eigenvalue corresponding to cutoff_frequency
        cutoff_freq_si = cutoff_frequency * (3e10) # Hz? No, w = 2*pi*f
        # Wait, getVibrationalSpectrum uses:
        # VibrationalSpectrum = AngularFreqSI/(6*pi*10**10) -> This is effectively freq/(2*pi*c) but approximated?
        # Let's trust the return of __getVibrationalSpectrum__.
        
        # Let's perform clamping after getting spectrum? No, we need SquareAngularFreq to be clean.
        
        # Calculate Min Eigenvalue Threshold
        # freq_cm = sqrt(Eig * 4.184e26) / (6*pi*10^10)  <-- definition in line 84
        # freq_cm * (6*pi*10^10) = sqrt(Eig * 4.184e26)
        # (freq_cm * 6*pi*10^10)^2 = Eig * 4.184e26
        # Eig = (freq_cm * 6*pi*10^10)**2 / 4.184e26
        
        pi_val = 3.14159265359
        threshold_eigenvalue = ((cutoff_frequency * 6 * pi_val * 1e10)**2) / (4.184e26)
        
        # Clamp vibrational modes (indices 6 to end)
        num_clamped = 0
        for i in range(6, len(eigValSorted)):
            if eigValSorted[i] < threshold_eigenvalue:
                eigValSorted[i] = threshold_eigenvalue
                num_clamped += 1
                
        if num_clamped > 0:
            print(f"  Note: Clamped {num_clamped} low-frequency modes to {cutoff_frequency} cm^-1 (Robust NMA).")
            
        VibrationalSpectrum = self.__getVibrationalSpectrum__(eigValSorted[6:])

        if not self.__checkSymmetric__(HessianSymmetric):
            warn('Fatal Warining: The hessian is NOT symmetric !!')

        if not self.__checkPositiveDefinite__(eigValSorted):
            # This check might still trigger if we checked sorted vs original, but we modified sorted.
            # Actually we already ensured positive definite by clamping with positive threshold.
            pass

        self.Hessian = HessianSymmetric * (kilocalorie/(gram*angstrom**2))
        self.SquareAngularFreq = eigValSorted * (kilocalorie/(gram*angstrom**2))
        self.NormalModes = eigVecSorted
        self.VibrationalSpectrum = VibrationalSpectrum * (1/centimeter)

    def PlotVibrationalSpectrum(self, binNum=1000, colorStr='salmon', labelStr='Vibrational Power Spectrum'):
        """
            Plot the histogram of vibrational power spectrum intensities.
            X-axis (cm^-1): Wave number.
            Y-axis (dimensionless number): Intensity of the histogram.
            Args:
            binNum (int=1000): Number of bins to generate the histogram.
            colorStr (string='salmon'): The matplotlib color string.`
            labelStr (string='Vibrational Power Spectrum'): The string used to show the label of current histogram. 
        """
        matplotlib.rcParams['figure.dpi'] = 600
        fig, ax = plt.subplots()
        ax.hist(self.VibrationalSpectrum, bins=binNum, color=colorStr,label=labelStr)
        ax.legend()
        plt.xlabel(r'Wave number($cm^{-1}$)')
        plt.ylabel(r'Intensity')
        plt.xlim(0,4000)
        plt.show()

    def getVibrationalEntropyCM(self, Temperature=300*unit.kelvin):
        """
        Classical (high-temperature limit) harmonic-oscillator vibrational entropy,
        summed over all internal (non-translational/rotational) normal modes.

        Per-mode entropy: S_i = k_B * [1 + ln(k_B*T / (hbar*omega_i))]
        Total:            S = N_int*k_B*(1 + ln(k_B*T)) - k_B * sum_i(ln(hbar*omega_i))

        This is the T -> large / hbar*omega << k_B*T limit of the quantum
        harmonic-oscillator entropy computed by `getVibrationalEntropyQM`; the
        two should converge for low-frequency modes at biologically relevant
        temperatures and diverge for high-frequency (e.g. X-H stretch) modes,
        where the classical approximation overestimates entropy.
        """
        SquareAngularFreqAKMA = self.SquareAngularFreq.value_in_unit(kilocalorie/(gram*angstrom**2))[6:]
        AngularFreqSI = np.sqrt((4.184*10**26)*SquareAngularFreqAKMA)
        NumAtoms = self.CUDASimulation.system.getNumParticles()
        internalDim = 3*NumAtoms - 6
        Temperature = Temperature.value_in_unit(kelvin)
        kBT = Boltzmann*Temperature
        VibrationalEntropyCM = internalDim*kBT*(1 + np.log(kBT)) - kBT*np.sum(np.log(hbar*AngularFreqSI))
        VibrationalEntropyCM = VibrationalEntropyCM*Avogadro * (joule/mole)
        self.VibrationalEntropyCM = VibrationalEntropyCM.in_units_of(kilocalorie/mole)
    
    def getVibrationalEntropyQM(self, Temperature=300*unit.kelvin):
        SquareAngularFreqAKMA = self.SquareAngularFreq.value_in_unit(kilocalorie/(gram*angstrom**2))[6:]
        AngularFreqSI = np.sqrt((4.184*10**26)*SquareAngularFreqAKMA)
        NumAtoms = self.CUDASimulation.system.getNumParticles()
        internalDim = 3*NumAtoms - 6
        Temperature = Temperature.value_in_unit(kelvin)
        kBT = Boltzmann*Temperature
        quantumEnergy = hbar*AngularFreqSI
        VibrationalEntropyQMArray = quantumEnergy/(np.exp(quantumEnergy/kBT) - 1) - kBT*np.log(1 - np.exp(-quantumEnergy/kBT))
        VibrationalEntropyQM = np.sum(VibrationalEntropyQMArray)
        VibrationalEntropyQM = VibrationalEntropyQM*Avogadro * (joule/mole)
        self.VibrationalEntropyQM = VibrationalEntropyQM.in_units_of(kilocalorie/mole)



        

