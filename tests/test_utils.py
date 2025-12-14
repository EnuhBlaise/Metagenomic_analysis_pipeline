import pytest
import tempfile
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

from metagenomics_toolkit.core.utils import (
    run_command, 
    ProcessError,
    create_conda_script,
    execute_conda_commands,
    ensure_directory,
    validate_file_exists,
    validate_directory_exists,
    find_files_with_extension,
    cleanup_temp_files
)

class TestRunCommand:
    """Test command execution utilities"""
    
    def test_successful_command(self):
        """Test successful command execution"""
        result = run_command("echo 'test'")
        assert result.returncode == 0
        assert "test" in result.stdout
    
    def test_failed_command(self):
        """Test failed command execution"""
        with pytest.raises(ProcessError):
            run_command("false")  # Command that always fails
    
    def test_command_with_args_list(self):
        """Test command execution with arguments list"""
        result = run_command(["echo", "test", "message"])
        assert result.returncode == 0
        assert "test message" in result.stdout
    
    def test_command_timeout(self):
        """Test command timeout"""
        with pytest.raises(ProcessError, match="timed out"):
            run_command("sleep 10", timeout=1)

class TestCondaScript:
    """Test conda script utilities"""
    
    def test_create_conda_script(self):
        """Test conda script creation"""
        script = create_conda_script(
            "/opt/conda/etc/profile.d/conda.sh",
            "test_env",
            ["echo 'hello'", "echo 'world'"]
        )
        
        assert "#!/bin/bash" in script
        assert "source /opt/conda/etc/profile.d/conda.sh" in script
        assert "conda activate test_env" in script
        assert "echo 'hello'" in script
        assert "echo 'world'" in script
        assert "conda deactivate" in script
    
    @patch('metagenomics_toolkit.core.utils.run_command')
    def test_execute_conda_commands(self, mock_run):
        """Test conda command execution"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_run.return_value = mock_result
        
        execute_conda_commands(
            "/opt/conda/etc/profile.d/conda.sh",
            "test_env",
            ["echo 'test'"]
        )
        
        # Verify run_command was called
        mock_run.assert_called_once()
        args, kwargs = mock_run.call_args
        assert args[0][0] == 'bash'  # First argument should be bash

class TestFileOperations:
    """Test file and directory utilities"""
    
    def test_ensure_directory(self):
        """Test directory creation"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_dir = Path(tmpdir) / "test" / "nested" / "directory"
            
            result = ensure_directory(test_dir)
            assert result.exists()
            assert result.is_dir()
            assert result == test_dir
    
    def test_validate_file_exists(self):
        """Test file validation"""
        with tempfile.NamedTemporaryFile() as tmp_file:
            # Should not raise exception
            result = validate_file_exists(tmp_file.name)
            assert result == Path(tmp_file.name)
        
        # Should raise exception for non-existent file
        with pytest.raises(FileNotFoundError):
            validate_file_exists("/nonexistent/file.txt")
    
    def test_validate_directory_exists(self):
        """Test directory validation"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Should not raise exception
            result = validate_directory_exists(tmpdir)
            assert result == Path(tmpdir)
        
        # Should raise exception for non-existent directory
        with pytest.raises(NotADirectoryError):
            validate_directory_exists("/nonexistent/directory")
    
    def test_find_files_with_extension(self):
        """Test finding files with specific extension"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            # Create test files
            (tmpdir_path / "test1.txt").touch()
            (tmpdir_path / "test2.txt").touch()
            (tmpdir_path / "test3.log").touch()
            
            # Create subdirectory with files
            subdir = tmpdir_path / "subdir"
            subdir.mkdir()
            (subdir / "test4.txt").touch()
            
            # Test non-recursive search
            txt_files = find_files_with_extension(tmpdir, "txt", recursive=False)
            assert len(txt_files) == 2
            assert all(f.suffix == ".txt" for f in txt_files)
            
            # Test recursive search
            txt_files_recursive = find_files_with_extension(tmpdir, "txt", recursive=True)
            assert len(txt_files_recursive) == 3
            
            # Test with dot prefix
            txt_files_dot = find_files_with_extension(tmpdir, ".txt", recursive=False)
            assert len(txt_files_dot) == 2
    
    def test_cleanup_temp_files(self):
        """Test temporary file cleanup"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            # Create test files and directories
            test_file = tmpdir_path / "test_file.txt"
            test_dir = tmpdir_path / "test_dir"
            test_file.touch()
            test_dir.mkdir()
            
            assert test_file.exists()
            assert test_dir.exists()
            
            # Clean up
            cleanup_temp_files(test_file, test_dir)
            
            assert not test_file.exists()
            assert not test_dir.exists()