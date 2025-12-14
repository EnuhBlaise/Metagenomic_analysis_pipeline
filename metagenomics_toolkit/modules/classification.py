import logging
from pathlib import Path
from typing import Optional

from ..core import utils, config

logger = logging.getLogger(__name__)

class GTDBTKClassifier:
    """GTDB-tk taxonomic classifier wrapper"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
        self.conda_script = self.config.get("conda.activate_script")
    
    def run(self, bins_dir: str, output_dir: str, 
            extension: str = "fa", **kwargs) -> str:
        """Run GTDB-tk classification
        
        Args:
            bins_dir: Directory containing genome bins
            output_dir: Output directory
            extension: File extension for genome files
            **kwargs: Additional parameters
            
        Returns:
            Path to classification results
        """
        logger.info(f"Starting GTDB-tk classification: {bins_dir} -> {output_dir}")
        
        # Validate inputs
        utils.validate_directory_exists(bins_dir, "Bins directory")
        utils.ensure_directory(output_dir)
        
        # Check if any genome files exist
        genome_files = utils.find_files_with_extension(bins_dir, extension)
        if not genome_files:
            raise FileNotFoundError(f"No genome files with extension '.{extension}' found in {bins_dir}")
        
        logger.info(f"Found {len(genome_files)} genome files for classification")
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("gtdbtk")
        threads = kwargs.get('threads', tool_config.get('default_threads', 8))
        
        conda_env = self.config.get_conda_env("gtdbtk")
        if not conda_env:
            raise ValueError("GTDB-tk conda environment not configured")
        
        # Build GTDB-tk command
        commands = [
            f"gtdbtk classify_wf -x {extension} --genome_dir {bins_dir} --out_dir {output_dir} --cpus {threads}"
        ]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Verify output files
            summary_files = [
                Path(output_dir) / "gtdbtk.bac120.summary.tsv",
                Path(output_dir) / "gtdbtk.ar53.summary.tsv"
            ]
            
            found_summaries = [f for f in summary_files if f.exists()]
            if not found_summaries:
                logger.warning("No GTDB-tk summary files found")
            else:
                logger.info(f"GTDB-tk classification completed. Found summary files: {[str(f) for f in found_summaries]}")
            
            return str(output_dir)
            
        except Exception as e:
            logger.error(f"GTDB-tk classification failed: {e}")
            raise

def run(args, config_manager: config.ConfigManager):
    """Main entry point for classification module"""
    
    if args.classify_tool == 'gtdbtk':
        classifier = GTDBTKClassifier(config_manager)
        results_dir = classifier.run(
            args.bins_dir,
            args.output,
            extension=getattr(args, 'extension', 'fa'),
            threads=getattr(args, 'threads', None)
        )
    else:
        raise ValueError(f"Unknown classification tool: {args.classify_tool}")
    
    logger.info(f"Classification completed successfully: {results_dir}")
    return results_dir