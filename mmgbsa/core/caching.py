"""
Caching system for MM/GBSA calculations.
Manages storage and retrieval of parameterized systems to accelerate repeated runs.
"""

import os
import json
import hashlib
import pickle
from pathlib import Path
import openmm as mm
import openmm.app as app


class CacheManager:
    """
    Manages caching of system data, topologies, and parameters.
    Accelerates repeated runs by storing compiled systems.
    """
    
    def __init__(self, cache_dir='.mmgbsa_cache', enabled=True, verbose=False):
        """
        Initialize cache manager.
        
        Parameters:
        -----------
        cache_dir : str
            Directory for cache storage
        enabled : bool
            Enable/disable caching
        verbose : bool
            Enable verbose logging
        """
        self.cache_dir = cache_dir
        self.enabled = enabled
        self.verbose = verbose
        
        if self.enabled:
            os.makedirs(self.cache_dir, exist_ok=True)
            if self.verbose:
                print(f"✓ Cache directory: {self.cache_dir}")
    
    def _get_file_hash(self, file_path):
        """
        Compute SHA256 hash of a file.
        
        Parameters:
        -----------
        file_path : str
            Path to the file
            
        Returns:
        --------
        str : Hex digest of file contents
        """
        sha256_hash = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    
    def get_cache_filename(self, topology_path, ligand_path=None, prefix='system'):
        """
        Generate cache filename based on input paths.
        
        Parameters:
        -----------
        topology_path : str
            Path to topology file
        ligand_path : str, optional
            Path to ligand file
        prefix : str
            Prefix for cache filename
            
        Returns:
        --------
        str : Full cache file path
        """
        if not self.enabled:
            return None
        
        try:
            # Hash topology
            topo_hash = self._get_file_hash(topology_path)[:8]
            
            # Hash ligand if provided
            if ligand_path and os.path.exists(ligand_path):
                ligand_hash = self._get_file_hash(ligand_path)[:8]
                cache_id = f"{prefix}_{topo_hash}_{ligand_hash}.pkl"
            else:
                cache_id = f"{prefix}_{topo_hash}.pkl"
            
            return os.path.join(self.cache_dir, cache_id)
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error computing cache filename: {e}")
            return None
    
    def save_system_to_cache(self, system, topology, cache_path):
        """
        Save OpenMM system and topology to cache.
        
        Parameters:
        -----------
        system : openmm.System
            OpenMM system object
        topology : openmm.app.Topology
            OpenMM topology object
        cache_path : str
            Path to save cache file
            
        Returns:
        --------
        bool : Success status
        """
        if not self.enabled or not cache_path:
            return False
        
        try:
            cache_data = {
                'system_xml': mm.XmlSerializer.serialize(system),
                'topology_xml': self._serialize_topology(topology)
            }
            
            with open(cache_path, 'wb') as f:
                pickle.dump(cache_data, f)
            
            if self.verbose:
                print(f"✓ Cached system to {cache_path}")
            return True
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error saving cache: {e}")
            return False
    
    def load_system_from_cache(self, cache_path):
        """
        Load OpenMM system and topology from cache.
        
        Parameters:
        -----------
        cache_path : str
            Path to cache file
            
        Returns:
        --------
        tuple : (system, topology) or (None, None) if not found/error
        """
        if not self.enabled or not cache_path or not os.path.exists(cache_path):
            return None, None
        
        try:
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)
            
            system = mm.XmlSerializer.deserialize(cache_data['system_xml'])
            topology = self._deserialize_topology(cache_data['topology_xml'])
            
            if self.verbose:
                print(f"✓ Loaded system from cache: {cache_path}")
            return system, topology
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error loading cache: {e}")
            return None, None
    
    def _serialize_topology(self, topology):
        """Serialize topology to XML string."""
        from io import StringIO
        with StringIO() as f:
            app.PDBFile.writeFile(topology, None, f)
            return f.getvalue()
    
    def _deserialize_topology(self, topology_xml):
        """Deserialize topology from XML string."""
        from io import StringIO
        pdb_file = app.PDBFile(StringIO(topology_xml))
        return pdb_file.topology
    
    def clear_cache(self):
        """Clear all cached files."""
        if not self.enabled:
            return
        
        try:
            import shutil
            if os.path.exists(self.cache_dir):
                shutil.rmtree(self.cache_dir)
                if self.verbose:
                    print(f"✓ Cleared cache directory: {self.cache_dir}")
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error clearing cache: {e}")
    
    def list_cache(self):
        """List all cached files."""
        if not self.enabled or not os.path.exists(self.cache_dir):
            return []
        
        try:
            return os.listdir(self.cache_dir)
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error listing cache: {e}")
            return []
