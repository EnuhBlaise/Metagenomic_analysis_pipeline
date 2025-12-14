import logging
from pathlib import Path
from typing import Optional

from ..core import utils, config

logger = logging.getLogger(__name__)

class CheckMQualityAssessor:
    """CheckM quality assessment tool wrapper"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
        self.conda_script = self.config.get("conda.activate_script")
    
    def run(self, bins_dir: str, output_dir: str, 
            extension: str = "fa", **kwargs) -> str:
        """Run CheckM quality assessment
        
        Args:
            bins_dir: Directory containing genome bins
            output_dir: Output directory
            extension: File extension for genome files
            **kwargs: Additional parameters
            
        Returns:
            Path to quality assessment results
        """
        logger.info(f"Starting CheckM quality assessment: {bins_dir} -> {output_dir}")
        
        # Validate inputs
        utils.validate_directory_exists(bins_dir, "Bins directory")
        utils.ensure_directory(output_dir)
        
        # Check if any genome files exist
        genome_files = utils.find_files_with_extension(bins_dir, extension)
        if not genome_files:
            raise FileNotFoundError(f"No genome files with extension '.{extension}' found in {bins_dir}")
        
        logger.info(f"Found {len(genome_files)} genome files for quality assessment")
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("checkm")
        threads = kwargs.get('threads', tool_config.get('default_threads', 8))
        
        conda_env = self.config.get_conda_env("checkm")
        if not conda_env:
            raise ValueError("CheckM conda environment not configured")
        
        # Build CheckM command
        commands = [
            f"checkm lineage_wf -t {threads} -x {extension} {bins_dir} {output_dir}"
        ]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Verify output files
            checkm_output_file = Path(output_dir) / "storage" / "bin_stats_ext.tsv"
            if not checkm_output_file.exists():
                # Try alternative output location
                checkm_output_file = Path(output_dir) / "lineage.ms"
            
            if checkm_output_file.exists():
                logger.info(f"CheckM quality assessment completed: {checkm_output_file}")
            else:
                logger.warning("CheckM output files not found in expected locations")
            
            return str(output_dir)
            
        except Exception as e:
            logger.error(f"CheckM quality assessment failed: {e}")
            raise

    def run_qa(self, checkm_dir: str, output_file: str = None, **kwargs) -> str:
        """Run CheckM qa command to generate summary
        
        Args:
            checkm_dir: CheckM output directory
            output_file: Optional output file for results
            **kwargs: Additional parameters
            
        Returns:
            Path to QA results
        """
        logger.info(f"Running CheckM qa: {checkm_dir}")
        
        utils.validate_directory_exists(checkm_dir, "CheckM directory")
        
        if output_file is None:
            output_file = Path(checkm_dir) / "checkm_qa_results.tsv"
        
        conda_env = self.config.get_conda_env("checkm")
        if not conda_env:
            raise ValueError("CheckM conda environment not configured")
        
        # Look for lineage.ms file
        lineage_file = Path(checkm_dir) / "lineage.ms"
        if not lineage_file.exists():
            raise FileNotFoundError(f"CheckM lineage file not found: {lineage_file}")
        
        # Build CheckM qa command
        commands = [
            f"checkm qa {lineage_file} {checkm_dir} --file {output_file} --tab_table"
        ]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            if Path(output_file).exists():
                logger.info(f"CheckM QA completed: {output_file}")
            else:
                logger.warning(f"CheckM QA output not found: {output_file}")
            
            return str(output_file)
            
        except Exception as e:
            logger.error(f"CheckM QA failed: {e}")
            raise

def run(args, config_manager: config.ConfigManager):
    """Main entry point for quality control module"""
    
    if args.qc_tool == 'checkm':
        assessor = CheckMQualityAssessor(config_manager)
        
        if hasattr(args, 'qa_mode') and args.qa_mode:
            # Run QA mode
            results_file = assessor.run_qa(
                args.bins_dir,
                output_file=getattr(args, 'output_file', None),
                threads=getattr(args, 'threads', None)
            )
        else:
            # Run full CheckM workflow
            results_dir = assessor.run(
                args.bins_dir,
                args.output,
                extension=getattr(args, 'extension', 'fa'),
                threads=getattr(args, 'threads', None)
            )
            results_file = results_dir
        
    else:
        raise ValueError(f"Unknown quality control tool: {args.qc_tool}")
    
    logger.info(f"Quality control completed successfully: {results_file}")
    return results_file