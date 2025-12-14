import subprocess
import os
import shutil
import signal
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
import logging

logger = logging.getLogger(__name__)

class ProcessError(Exception):
    """Custom exception for process execution errors"""
    pass

def run_command(cmd: Union[str, List[str]], 
                cwd: Optional[str] = None,
                env: Optional[Dict[str, str]] = None,
                timeout: Optional[int] = None,
                capture_output: bool = True) -> subprocess.CompletedProcess:
    """Run shell command with error handling
    
    Args:
        cmd: Command to run (string or list)
        cwd: Working directory
        env: Environment variables
        timeout: Timeout in seconds
        capture_output: Whether to capture stdout/stderr
        
    Returns:
        CompletedProcess instance
        
    Raises:
        ProcessError: If command fails
    """
    if isinstance(cmd, str):
        cmd = cmd.split()
    
    logger.info(f"Running command: {' '.join(cmd)}")
    
    try:
        # Merge environment variables
        full_env = os.environ.copy()
        if env:
            full_env.update(env)
        
        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=full_env,
            timeout=timeout,
            capture_output=capture_output,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            error_msg = f"Command failed with return code {result.returncode}"
            if capture_output and result.stderr:
                error_msg += f"\nSTDERR: {result.stderr}"
            logger.error(error_msg)
            raise ProcessError(error_msg)
        
        if capture_output and result.stdout:
            logger.debug(f"STDOUT: {result.stdout}")
            
        return result
        
    except subprocess.TimeoutExpired:
        error_msg = f"Command timed out after {timeout} seconds"
        logger.error(error_msg)
        raise ProcessError(error_msg)
    except Exception as e:
        error_msg = f"Error running command: {e}"
        logger.error(error_msg)
        raise ProcessError(error_msg)

def create_conda_script(conda_activate_script: str, 
                       environment: str, 
                       commands: List[str]) -> str:
    """Create bash script with conda environment activation
    
    Args:
        conda_activate_script: Path to conda activation script
        environment: Conda environment name or path
        commands: List of commands to run in environment
        
    Returns:
        Bash script content
    """
    script_lines = [
        "#!/bin/bash",
        "set -e",
        "",
        f"source {conda_activate_script}",
        "unset PYTHONPATH",
        "",
        f"conda activate {environment}",
        ""
    ]
    
    script_lines.extend(commands)
    script_lines.extend([
        "",
        "conda deactivate"
    ])
    
    return "\n".join(script_lines)

def execute_conda_commands(conda_activate_script: str,
                          environment: str,
                          commands: List[str],
                          cwd: Optional[str] = None,
                          timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """Execute commands in a conda environment
    
    Args:
        conda_activate_script: Path to conda activation script
        environment: Conda environment name or path
        commands: List of commands to run
        cwd: Working directory
        timeout: Timeout in seconds
        
    Returns:
        CompletedProcess instance
    """
    script_content = create_conda_script(conda_activate_script, environment, commands)
    
    # Create temporary script file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
        f.write(script_content)
        script_path = f.name
    
    try:
        # Make script executable
        os.chmod(script_path, 0o755)
        
        # Execute script
        result = run_command(['bash', script_path], cwd=cwd, timeout=timeout)
        return result
        
    finally:
        # Clean up temporary script
        try:
            os.unlink(script_path)
        except OSError:
            pass

def ensure_directory(path: Union[str, Path]) -> Path:
    """Ensure directory exists, create if necessary
    
    Args:
        path: Directory path
        
    Returns:
        Path object
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def validate_file_exists(file_path: Union[str, Path], description: str = "File") -> Path:
    """Validate that a file exists
    
    Args:
        file_path: Path to file
        description: Description for error messages
        
    Returns:
        Path object
        
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")
    if not path.is_file():
        raise ValueError(f"{description} is not a file: {path}")
    return path

def validate_directory_exists(dir_path: Union[str, Path], description: str = "Directory") -> Path:
    """Validate that a directory exists
    
    Args:
        dir_path: Path to directory
        description: Description for error messages
        
    Returns:
        Path object
        
    Raises:
        NotADirectoryError: If directory doesn't exist
    """
    path = Path(dir_path)
    if not path.exists():
        raise NotADirectoryError(f"{description} not found: {path}")
    if not path.is_dir():
        raise ValueError(f"{description} is not a directory: {path}")
    return path

def find_files_with_extension(directory: Union[str, Path], 
                             extension: str,
                             recursive: bool = False) -> List[Path]:
    """Find files with specific extension in directory
    
    Args:
        directory: Directory to search
        extension: File extension (with or without dot)
        recursive: Whether to search recursively
        
    Returns:
        List of Path objects
    """
    directory = Path(directory)
    
    # Ensure extension starts with dot
    if not extension.startswith('.'):
        extension = '.' + extension
    
    pattern = f"*{extension}"
    if recursive:
        return list(directory.rglob(pattern))
    else:
        return list(directory.glob(pattern))

def get_file_size_mb(file_path: Union[str, Path]) -> float:
    """Get file size in megabytes
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in MB
    """
    path = Path(file_path)
    return path.stat().st_size / (1024 * 1024)

def check_disk_space(directory: Union[str, Path], required_gb: float = 10.0) -> bool:
    """Check if sufficient disk space is available
    
    Args:
        directory: Directory to check
        required_gb: Required space in GB
        
    Returns:
        True if sufficient space available
    """
    try:
        statvfs = shutil.disk_usage(directory)
        available_gb = statvfs.free / (1024 ** 3)
        return available_gb >= required_gb
    except Exception:
        return False

def cleanup_temp_files(*file_paths: Union[str, Path]):
    """Clean up temporary files safely
    
    Args:
        *file_paths: Variable number of file paths to clean up
    """
    for file_path in file_paths:
        try:
            path = Path(file_path)
            if path.exists():
                if path.is_file():
                    path.unlink()
                elif path.is_dir():
                    shutil.rmtree(path)
                logger.debug(f"Cleaned up: {path}")
        except Exception as e:
            logger.warning(f"Could not clean up {file_path}: {e}")

class TimeoutHandler:
    """Context manager for handling timeouts with signal"""
    
    def __init__(self, timeout: int):
        self.timeout = timeout
        self.old_handler = None
    
    def __enter__(self):
        if self.timeout:
            self.old_handler = signal.signal(signal.SIGALRM, self._timeout_handler)
            signal.alarm(self.timeout)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.timeout:
            signal.alarm(0)
            if self.old_handler:
                signal.signal(signal.SIGALRM, self.old_handler)
    
    def _timeout_handler(self, signum, frame):
        raise TimeoutError(f"Operation timed out after {self.timeout} seconds")