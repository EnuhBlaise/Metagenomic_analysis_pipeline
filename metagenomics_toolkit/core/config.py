import json
import os
from pathlib import Path
from typing import Dict, Any, Optional

class ConfigManager:
    """Manages configuration for the metagenomics toolkit"""
    
    DEFAULT_CONFIG = {
        "conda": {
            "activate_script": "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh",
            "environments": {
                "flye": "/opt/bifxapps/miniconda3/envs/flye",
                "metabat2": "/home/glbrc.org/benuh/.conda/envs/metabat2_env",
                "gtdbtk": "/home/glbrc.org/benuh/.conda/envs/gtdbtk_env",
                "checkm": "/opt/bifxapps/miniconda3/envs/checkm",
                "prokka": "/opt/bifxapps/miniconda3/envs/prokka-1.14.6",
                "metawrap": "/opt/bifxapps/miniconda3/envs/metawrap",
                "racon": "/home/glbrc.org/benuh/.conda/envs/racon_env",
                "minimap2": "/home/glbrc.org/benuh/.conda/envs/minimap2_env",
                "samtools": "/home/glbrc.org/benuh/.conda/envs/samtools_env"
            }
        },
        "tools": {
            "flye": {
                "default_threads": 20,
                "meta_mode": True,
                "plasmids": True,
                "debug": True
            },
            "metabat2": {
                "min_contig_length": 2000,
                "default_threads": 8
            },
            "gtdbtk": {
                "default_threads": 8,
                "file_extension": "fa"
            },
            "checkm": {
                "default_threads": 8,
                "file_extension": "fa"
            },
            "prokka": {
                "force_overwrite": True
            },
            "spades": {
                "meta_mode": True,
                "default_threads": 4
            },
            "megahit": {
                "meta_mode": True,
                "default_threads": 4
            }
        },
        "paths": {
            "temp_dir": "/tmp/metagenomics_toolkit",
            "script_dir": None  # Will be set dynamically
        },
        "logging": {
            "level": "INFO",
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        }
    }
    
    def __init__(self, config_file: Optional[str] = None):
        """Initialize configuration manager
        
        Args:
            config_file: Path to custom configuration file
        """
        self.config = self.DEFAULT_CONFIG.copy()
        self._set_script_directory()
        
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
            
    def _set_script_directory(self):
        """Set the script directory to where original bash scripts are located"""
        current_dir = Path(__file__).parent.parent.parent
        self.config["paths"]["script_dir"] = str(current_dir)
    
    def load_config(self, config_file: str):
        """Load configuration from JSON file"""
        try:
            with open(config_file, 'r') as f:
                custom_config = json.load(f)
            self._merge_config(custom_config)
        except Exception as e:
            raise ValueError(f"Error loading config file {config_file}: {e}")
    
    def _merge_config(self, custom_config: Dict[str, Any]):
        """Recursively merge custom configuration with default"""
        def merge_dict(base: Dict, update: Dict):
            for key, value in update.items():
                if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                    merge_dict(base[key], value)
                else:
                    base[key] = value
        
        merge_dict(self.config, custom_config)
    
    def get(self, key_path: str, default=None):
        """Get configuration value using dot notation
        
        Args:
            key_path: Dot-separated path to config value (e.g., 'conda.environments.flye')
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        keys = key_path.split('.')
        value = self.config
        
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
                
        return value
    
    def get_conda_env(self, tool: str) -> str:
        """Get conda environment for specific tool"""
        return self.get(f"conda.environments.{tool}")
    
    def get_tool_config(self, tool: str) -> Dict[str, Any]:
        """Get tool-specific configuration"""
        return self.get(f"tools.{tool}", {})
    
    def save_config(self, config_file: str):
        """Save current configuration to file"""
        try:
            with open(config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
        except Exception as e:
            raise ValueError(f"Error saving config file {config_file}: {e}")
    
    def create_default_config(self, config_file: str):
        """Create a default configuration file"""
        self.save_config(config_file)