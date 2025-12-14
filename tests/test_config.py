import pytest
import json
import tempfile
from pathlib import Path

from metagenomics_toolkit.core.config import ConfigManager

class TestConfigManager:
    """Test configuration management"""
    
    def test_default_config(self):
        """Test default configuration loading"""
        config = ConfigManager()
        
        # Check that default config is loaded
        assert config.get("conda.activate_script") is not None
        assert config.get("tools.flye.default_threads") == 20
        assert config.get("tools.metabat2.min_contig_length") == 2000
    
    def test_get_method(self):
        """Test config get method with dot notation"""
        config = ConfigManager()
        
        # Test existing keys
        assert config.get("conda.activate_script") == "/opt/bifxapps/miniconda3/etc/profile.d/conda.sh"
        assert config.get("tools.flye.meta_mode") == True
        
        # Test default values
        assert config.get("nonexistent.key", "default") == "default"
        assert config.get("another.nonexistent.key") is None
    
    def test_conda_env_method(self):
        """Test getting conda environments"""
        config = ConfigManager()
        
        flye_env = config.get_conda_env("flye")
        assert flye_env == "/opt/bifxapps/miniconda3/envs/flye"
        
        nonexistent_env = config.get_conda_env("nonexistent")
        assert nonexistent_env is None
    
    def test_tool_config_method(self):
        """Test getting tool configurations"""
        config = ConfigManager()
        
        flye_config = config.get_tool_config("flye")
        assert isinstance(flye_config, dict)
        assert flye_config["default_threads"] == 20
        
        nonexistent_config = config.get_tool_config("nonexistent")
        assert nonexistent_config == {}
    
    def test_custom_config_loading(self):
        """Test loading custom configuration file"""
        custom_config = {
            "conda": {
                "environments": {
                    "flye": "/custom/path/to/flye"
                }
            },
            "tools": {
                "flye": {
                    "default_threads": 50
                }
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(custom_config, f)
            config_file = f.name
        
        try:
            config = ConfigManager(config_file)
            
            # Check that custom config overrides defaults
            assert config.get_conda_env("flye") == "/custom/path/to/flye"
            assert config.get("tools.flye.default_threads") == 50
            
            # Check that non-overridden defaults remain
            assert config.get("tools.flye.meta_mode") == True
            
        finally:
            Path(config_file).unlink()
    
    def test_config_save_and_load(self):
        """Test saving and loading configuration"""
        config = ConfigManager()
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            config_file = f.name
        
        try:
            # Save config
            config.save_config(config_file)
            
            # Load it again
            new_config = ConfigManager(config_file)
            
            # Should be identical
            assert new_config.config == config.config
            
        finally:
            Path(config_file).unlink()
    
    def test_create_default_config(self):
        """Test creating default configuration file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            config_file = f.name
        
        try:
            config = ConfigManager()
            config.create_default_config(config_file)
            
            # Verify file was created and is valid JSON
            assert Path(config_file).exists()
            
            with open(config_file) as f:
                loaded_config = json.load(f)
            
            assert isinstance(loaded_config, dict)
            assert "conda" in loaded_config
            assert "tools" in loaded_config
            
        finally:
            Path(config_file).unlink()