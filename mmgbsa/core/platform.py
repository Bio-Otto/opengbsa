"""
Platform configuration and management for MM/GBSA calculations.
Handles OpenMM platform selection (CPU, CUDA, OpenCL) for main analysis and decomposition.
"""

import openmm as mm


class PlatformManager:
    """
    Manages OpenMM platform selection and configuration.
    Supports separate platforms for main analysis vs decomposition.
    """
    
    def __init__(self, verbose=False):
        """
        Initialize platform manager.
        
        Parameters:
        -----------
        verbose : bool
            Enable verbose logging
        """
        self.verbose = verbose
        self.preferred_platform = None
        self.decomposition_platform = None
        self.prefer_cuda = False
        self.prefer_opencl = False
        
    def set_platform_settings(self, platform_settings):
        """
        Configure platform preferences for analysis and decomposition.
        
        Parameters:
        -----------
        platform_settings : dict
            Configuration dictionary with optional keys:
            - 'preferred_platform': Platform for main analysis (default: auto-detect)
            - 'decomposition_platform': Platform for decomposition (default: same as analysis)
            Example: {'preferred_platform': 'CUDA', 'decomposition_platform': 'CPU'}
        """
        if not platform_settings:
            return
        
        if isinstance(platform_settings, dict):
            self.preferred_platform = platform_settings.get('preferred_platform')
            self.decomposition_platform = platform_settings.get('decomposition_platform')
            
            if self.verbose:
                if self.preferred_platform:
                    print(f"  Platform for analysis: {self.preferred_platform}")
                if self.decomposition_platform:
                    print(f"  Platform for decomposition: {self.decomposition_platform}")
    
    def setup_optimized_platform(self):
        """
        Set up optimized OpenMM platform based on preferences and availability.
        
        Returns:
        --------
        openmm.Platform : Configured OpenMM platform or None
        """
        platform_name = None
        
        try:
            # Check available platforms
            available_platforms = [mm.Platform.getPlatform(i).getName() 
                                  for i in range(mm.Platform.getNumPlatforms())]
            
            if self.verbose:
                print(f"Available platforms: {available_platforms}")
            
            # Use preferred platform if specified and available
            if self.preferred_platform:
                if self.preferred_platform in available_platforms:
                    platform_name = self.preferred_platform
                else:
                    if self.verbose:
                        print(f"⚠ Requested platform '{self.preferred_platform}' not available. "
                              f"Available: {available_platforms}")
            
            # Auto-select if no preference
            if not platform_name:
                # Try CUDA first, then OpenCL, fall back to CPU
                for platform in ['CUDA', 'OpenCL', 'CPU']:
                    if platform in available_platforms:
                        platform_name = platform
                        break
            
            if platform_name:
                platform = mm.Platform.getPlatformByName(platform_name)
                if self.verbose:
                    print(f"Using platform: {platform_name}")
                return platform
            else:
                if self.verbose:
                    print("No suitable platform found")
                return None
                
        except Exception as e:
            if self.verbose:
                print(f"Error setting up platform: {e}")
            return None
    
    def get_preferred_platform(self):
        """Get the user's preferred platform name (if set)."""
        return self.preferred_platform
    
    def get_decomposition_platform(self):
        """Get the decomposition platform name (if set)."""
        return self.decomposition_platform
