import logging
from pathlib import Path
from typing import Optional, Dict, Any

from ..core import utils, config

logger = logging.getLogger(__name__)

class AssemblyRunner:
    """Base class for assembly tools"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
        self.conda_script = self.config.get("conda.activate_script")
    
    def validate_inputs(self, input_file: str, output_dir: str):
        """Validate input files and output directory"""
        utils.validate_file_exists(input_file, "Input FASTQ file")
        utils.ensure_directory(output_dir)

class FlyeAssembler(AssemblyRunner):
    """FLYE metagenomic assembler wrapper"""
    
    def run(self, input_file: str, output_dir: str, 
            read_type: str = "nano-raw", **kwargs) -> str:
        """Run FLYE assembly
        
        Args:
            input_file: Input FASTQ file
            output_dir: Output directory
            read_type: Type of reads (nano-raw, nano-corr, pacbio)
            **kwargs: Additional parameters
            
        Returns:
            Path to assembly file
        """
        logger.info(f"Starting FLYE assembly: {input_file} -> {output_dir}")
        
        self.validate_inputs(input_file, output_dir)
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("flye")
        threads = kwargs.get('threads', tool_config.get('default_threads', 20))
        
        # Build FLYE command
        cmd_parts = [
            "flye",
            f"--{read_type}", input_file,
            "-t", str(threads),
            "-o", output_dir
        ]
        
        # Add optional flags
        if tool_config.get('meta_mode', True):
            cmd_parts.append("--meta")
        if tool_config.get('plasmids', True):
            cmd_parts.append("--plasmids")
        if tool_config.get('debug', True):
            cmd_parts.append("--debug")
        
        # Execute command
        conda_env = self.config.get_conda_env("flye")
        if not conda_env:
            raise ValueError("FLYE conda environment not configured")
        
        commands = [" ".join(cmd_parts)]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Return path to assembly file
            assembly_file = Path(output_dir) / "assembly.fasta"
            if not assembly_file.exists():
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")
            
            logger.info(f"FLYE assembly completed: {assembly_file}")
            return str(assembly_file)
            
        except Exception as e:
            logger.error(f"FLYE assembly failed: {e}")
            raise

class SpadesAssembler(AssemblyRunner):
    """SPAdes metagenomic assembler wrapper"""
    
    def run(self, input_file: str, output_dir: str, **kwargs) -> str:
        """Run SPAdes assembly
        
        Args:
            input_file: Input FASTQ file(s) - can be comma-separated for paired-end
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Path to assembly file
        """
        logger.info(f"Starting SPAdes assembly: {input_file} -> {output_dir}")
        
        # Handle paired-end reads
        input_files = input_file.split(',')
        if len(input_files) == 1:
            # Single-end
            for file in input_files:
                utils.validate_file_exists(file.strip(), "Input FASTQ file")
        elif len(input_files) == 2:
            # Paired-end
            read1, read2 = [f.strip() for f in input_files]
            utils.validate_file_exists(read1, "Read 1 FASTQ file")
            utils.validate_file_exists(read2, "Read 2 FASTQ file")
        else:
            raise ValueError("SPAdes supports single-end or paired-end reads only")
        
        utils.ensure_directory(output_dir)
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("spades")
        threads = kwargs.get('threads', tool_config.get('default_threads', 4))
        
        # Build SPAdes command
        cmd_parts = [
            "spades.py",
            "-t", str(threads),
            "-o", output_dir
        ]
        
        if tool_config.get('meta_mode', True):
            cmd_parts.append("--meta")
        
        # Add read files
        if len(input_files) == 1:
            cmd_parts.extend(["-s", input_files[0].strip()])
        else:
            cmd_parts.extend(["-1", read1, "-2", read2])
        
        # Execute command
        conda_env = self.config.get_conda_env("spades")
        if not conda_env:
            # Try using system SPAdes or default environment
            conda_env = "base"
        
        commands = [" ".join(cmd_parts)]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Return path to assembly file
            assembly_file = Path(output_dir) / "contigs.fasta"
            if not assembly_file.exists():
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")
            
            logger.info(f"SPAdes assembly completed: {assembly_file}")
            return str(assembly_file)
            
        except Exception as e:
            logger.error(f"SPAdes assembly failed: {e}")
            raise

class MegahitAssembler(AssemblyRunner):
    """MEGAHIT metagenomic assembler wrapper"""
    
    def run(self, input_file: str, output_dir: str, **kwargs) -> str:
        """Run MEGAHIT assembly
        
        Args:
            input_file: Input FASTQ file(s)
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Path to assembly file
        """
        logger.info(f"Starting MEGAHIT assembly: {input_file} -> {output_dir}")
        
        # Handle paired-end reads
        input_files = input_file.split(',')
        if len(input_files) == 1:
            utils.validate_file_exists(input_files[0].strip(), "Input FASTQ file")
        elif len(input_files) == 2:
            read1, read2 = [f.strip() for f in input_files]
            utils.validate_file_exists(read1, "Read 1 FASTQ file")
            utils.validate_file_exists(read2, "Read 2 FASTQ file")
        else:
            raise ValueError("MEGAHIT supports single-end or paired-end reads only")
        
        utils.ensure_directory(output_dir)
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("megahit")
        threads = kwargs.get('threads', tool_config.get('default_threads', 4))
        
        # Build MEGAHIT command
        cmd_parts = [
            "megahit",
            "-t", str(threads),
            "-o", output_dir
        ]
        
        if tool_config.get('meta_mode', True):
            cmd_parts.append("--meta")
        
        # Add read files
        if len(input_files) == 1:
            cmd_parts.extend(["-r", input_files[0].strip()])
        else:
            cmd_parts.extend(["-1", read1, "-2", read2])
        
        # Execute command
        conda_env = self.config.get_conda_env("megahit")
        if not conda_env:
            conda_env = "base"
        
        commands = [" ".join(cmd_parts)]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Return path to assembly file
            assembly_file = Path(output_dir) / "final.contigs.fa"
            if not assembly_file.exists():
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")
            
            logger.info(f"MEGAHIT assembly completed: {assembly_file}")
            return str(assembly_file)
            
        except Exception as e:
            logger.error(f"MEGAHIT assembly failed: {e}")
            raise

def run(args, config_manager: config.ConfigManager):
    """Main entry point for assembly module"""
    
    if args.assembly_tool == 'flye':
        assembler = FlyeAssembler(config_manager)
        assembly_file = assembler.run(
            args.input, 
            args.output,
            read_type=getattr(args, 'read_type', 'nano-raw'),
            threads=getattr(args, 'threads', None)
        )
    elif args.assembly_tool == 'spades':
        assembler = SpadesAssembler(config_manager)
        assembly_file = assembler.run(
            args.input,
            args.output,
            threads=getattr(args, 'threads', None)
        )
    elif args.assembly_tool == 'megahit':
        assembler = MegahitAssembler(config_manager)
        assembly_file = assembler.run(
            args.input,
            args.output,
            threads=getattr(args, 'threads', None)
        )
    else:
        raise ValueError(f"Unknown assembly tool: {args.assembly_tool}")
    
    logger.info(f"Assembly completed successfully: {assembly_file}")
    return assembly_file