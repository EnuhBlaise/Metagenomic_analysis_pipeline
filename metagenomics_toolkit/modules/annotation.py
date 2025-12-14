import logging
from pathlib import Path
from typing import Optional, List

from ..core import utils, config

logger = logging.getLogger(__name__)

class ProkkaAnnotator:
    """Prokka genome annotation tool wrapper"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
        self.conda_script = self.config.get("conda.activate_script")
    
    def run(self, input_path: str, output_dir: str, **kwargs) -> str:
        """Run Prokka annotation
        
        Args:
            input_path: Input genome file or directory containing genomes
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Path to annotation results
        """
        logger.info(f"Starting Prokka annotation: {input_path} -> {output_dir}")
        
        # Determine if input is file or directory
        input_path_obj = Path(input_path)
        
        if input_path_obj.is_file():
            # Single genome file
            utils.validate_file_exists(input_path, "Input genome file")
            genome_files = [input_path_obj]
        elif input_path_obj.is_dir():
            # Directory with multiple genomes
            utils.validate_directory_exists(input_path, "Input genomes directory")
            # Find FASTA files
            extensions = ['fa', 'fasta', 'fna']
            genome_files = []
            for ext in extensions:
                genome_files.extend(utils.find_files_with_extension(input_path, ext))
            
            if not genome_files:
                raise FileNotFoundError(f"No genome files found in {input_path}")
        else:
            raise ValueError(f"Input path does not exist: {input_path}")
        
        utils.ensure_directory(output_dir)
        
        logger.info(f"Found {len(genome_files)} genome files for annotation")
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("prokka")
        force_overwrite = tool_config.get('force_overwrite', True)
        
        conda_env = self.config.get_conda_env("prokka")
        if not conda_env:
            raise ValueError("Prokka conda environment not configured")
        
        results = []
        
        # Process each genome file
        for i, genome_file in enumerate(genome_files):
            genome_name = genome_file.stem
            genome_output_dir = Path(output_dir) / f"prokka_{genome_name}"
            
            # Build Prokka command
            cmd_parts = [
                "prokka",
                "--outdir", str(genome_output_dir),
                str(genome_file)
            ]
            
            if force_overwrite:
                cmd_parts.append("--force")
            
            # Add optional parameters
            if 'prefix' in kwargs:
                cmd_parts.extend(["--prefix", kwargs['prefix']])
            else:
                cmd_parts.extend(["--prefix", genome_name])
            
            if 'genus' in kwargs:
                cmd_parts.extend(["--genus", kwargs['genus']])
            
            if 'species' in kwargs:
                cmd_parts.extend(["--species", kwargs['species']])
            
            commands = [" ".join(cmd_parts)]
            
            try:
                logger.info(f"Annotating genome {i+1}/{len(genome_files)}: {genome_file.name}")
                
                utils.execute_conda_commands(
                    self.conda_script,
                    conda_env,
                    commands,
                    timeout=kwargs.get('timeout')
                )
                
                # Verify output
                gff_file = genome_output_dir / f"{genome_name}.gff"
                if not gff_file.exists():
                    logger.warning(f"Prokka output not found for {genome_file.name}")
                else:
                    results.append(str(genome_output_dir))
                    logger.info(f"Prokka annotation completed for {genome_file.name}")
                
            except Exception as e:
                logger.error(f"Prokka annotation failed for {genome_file.name}: {e}")
                if kwargs.get('continue_on_error', False):
                    continue
                else:
                    raise
        
        if not results:
            raise RuntimeError("No successful Prokka annotations completed")
        
        logger.info(f"Prokka annotation completed for {len(results)} genomes")
        return str(output_dir)

def run(args, config_manager: config.ConfigManager):
    """Main entry point for annotation module"""
    
    if args.annotate_tool == 'prokka':
        annotator = ProkkaAnnotator(config_manager)
        results_dir = annotator.run(
            args.input,
            args.output,
            threads=getattr(args, 'threads', None),
            genus=getattr(args, 'genus', None),
            species=getattr(args, 'species', None),
            prefix=getattr(args, 'prefix', None),
            continue_on_error=getattr(args, 'continue_on_error', False)
        )
    else:
        raise ValueError(f"Unknown annotation tool: {args.annotate_tool}")
    
    logger.info(f"Annotation completed successfully: {results_dir}")
    return results_dir