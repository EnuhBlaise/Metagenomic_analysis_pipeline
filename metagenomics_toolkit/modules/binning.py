import logging
from pathlib import Path
from typing import Optional, List

from ..core import utils, config

logger = logging.getLogger(__name__)

class BinningRunner:
    """Base class for binning tools"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
        self.conda_script = self.config.get("conda.activate_script")
    
    def validate_inputs(self, contigs_file: str, reads_file: str, output_dir: str):
        """Validate input files and output directory"""
        utils.validate_file_exists(contigs_file, "Contigs file")
        utils.validate_file_exists(reads_file, "Reads file")
        utils.ensure_directory(output_dir)

class MetaBat2Binner(BinningRunner):
    """MetaBAT2 binning tool wrapper"""
    
    def run(self, contigs_file: str, reads_file: str, output_dir: str, **kwargs) -> str:
        """Run MetaBAT2 binning
        
        Args:
            contigs_file: Input contigs FASTA file
            reads_file: Input reads FASTQ file
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Path to binning directory
        """
        logger.info(f"Starting MetaBAT2 binning: {contigs_file} -> {output_dir}")
        
        self.validate_inputs(contigs_file, reads_file, output_dir)
        
        # Get tool configuration
        tool_config = self.config.get_tool_config("metabat2")
        min_contig_length = kwargs.get('min_length', tool_config.get('min_contig_length', 2000))
        threads = kwargs.get('threads', tool_config.get('default_threads', 8))
        
        # Create mapping and binning directories
        output_path = Path(output_dir)
        mapping_dir = output_path / "mapping_output"
        binning_dir = output_path / "binning_output"
        utils.ensure_directory(mapping_dir)
        utils.ensure_directory(binning_dir)
        
        conda_env = self.config.get_conda_env("metabat2")
        if not conda_env:
            raise ValueError("MetaBAT2 conda environment not configured")
        
        commands = []
        
        # Step 1: Map reads to contigs with minimap2
        minimap2_env = self.config.get_conda_env("minimap2")
        if minimap2_env:
            # Use separate minimap2 environment
            minimap2_script = utils.create_conda_script(
                self.conda_script,
                minimap2_env,
                [f"minimap2 -t {threads} -ax map-hifi {contigs_file} {reads_file} > {mapping_dir}/minimap2.sam"]
            )
            
            # Create and execute minimap2 script
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
                f.write(minimap2_script)
                minimap2_script_path = f.name
            
            try:
                utils.run_command(['bash', minimap2_script_path])
            finally:
                Path(minimap2_script_path).unlink()
        else:
            # Use MetaBAT2 environment for mapping
            commands.append(f"minimap2 -t {threads} -ax map-hifi {contigs_file} {reads_file} > {mapping_dir}/minimap2.sam")
        
        # Step 2: Convert SAM to sorted BAM
        samtools_env = self.config.get_conda_env("samtools")
        if samtools_env:
            samtools_script = utils.create_conda_script(
                self.conda_script,
                samtools_env,
                [
                    f"samtools view -bS {mapping_dir}/minimap2.sam > {mapping_dir}/minimap2.bam",
                    f"samtools sort -o {mapping_dir}/minimap2_sorted.bam {mapping_dir}/minimap2.bam"
                ]
            )
            
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
                f.write(samtools_script)
                samtools_script_path = f.name
            
            try:
                utils.run_command(['bash', samtools_script_path])
            finally:
                Path(samtools_script_path).unlink()
        else:
            commands.extend([
                f"samtools view -bS {mapping_dir}/minimap2.sam > {mapping_dir}/minimap2.bam",
                f"samtools sort -o {mapping_dir}/minimap2_sorted.bam {mapping_dir}/minimap2.bam"
            ])
        
        # Step 3: Run MetaBAT2
        commands.append(f"runMetaBat.sh -m {min_contig_length} {contigs_file} {mapping_dir}/minimap2_sorted.bam")
        
        # Change to binning output directory for MetaBAT2 output
        commands.insert(-1, f"cd {binning_dir}")
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                cwd=str(binning_dir),
                timeout=kwargs.get('timeout')
            )
            
            # Find MetaBAT2 output directory
            metabat_dirs = list(binning_dir.glob("*.metabat-bins*"))
            if not metabat_dirs:
                raise FileNotFoundError("MetaBAT2 output directory not found")
            
            bins_dir = metabat_dirs[0]
            logger.info(f"MetaBAT2 binning completed: {bins_dir}")
            return str(bins_dir)
            
        except Exception as e:
            logger.error(f"MetaBAT2 binning failed: {e}")
            raise

class ConcoTBinner(BinningRunner):
    """CONCOCT binning tool wrapper"""
    
    def run(self, contigs_file: str, reads_file: str, output_dir: str, **kwargs) -> str:
        """Run CONCOCT binning using MetaWRAP
        
        Args:
            contigs_file: Input contigs FASTA file
            reads_file: Input reads FASTQ file
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Path to binning directory
        """
        logger.info(f"Starting CONCOCT binning: {contigs_file} -> {output_dir}")
        
        self.validate_inputs(contigs_file, reads_file, output_dir)
        
        threads = kwargs.get('threads', 96)  # CONCOCT typically uses more threads
        
        conda_env = self.config.get_conda_env("metawrap")
        if not conda_env:
            raise ValueError("MetaWRAP conda environment not configured")
        
        # Build MetaWRAP binning command
        commands = [
            f"metawrap binning -o {output_dir} -t {threads} -a {contigs_file} --concoct {reads_file} --single-end --run-checkm"
        ]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            # Find CONCOCT output directory
            concoct_dir = Path(output_dir) / "CONCOCT_BINNING" / "concoct_bins"
            if not concoct_dir.exists():
                raise FileNotFoundError(f"CONCOCT output directory not found: {concoct_dir}")
            
            logger.info(f"CONCOCT binning completed: {concoct_dir}")
            return str(concoct_dir)
            
        except Exception as e:
            logger.error(f"CONCOCT binning failed: {e}")
            raise

class MetaWRAPBinner(BinningRunner):
    """MetaWRAP comprehensive binning wrapper"""
    
    def run(self, contigs_file: str, reads_file: str, output_dir: str, 
            binners: List[str] = None, **kwargs) -> str:
        """Run MetaWRAP binning with multiple tools
        
        Args:
            contigs_file: Input contigs FASTA file
            reads_file: Input reads FASTQ file
            output_dir: Output directory
            binners: List of binners to use (metabat1, metabat2, maxbin2, concoct)
            **kwargs: Additional parameters
            
        Returns:
            Path to binning directory
        """
        logger.info(f"Starting MetaWRAP binning: {contigs_file} -> {output_dir}")
        
        self.validate_inputs(contigs_file, reads_file, output_dir)
        
        if binners is None:
            binners = ['metabat2', 'maxbin2', 'concoct']
        
        threads = kwargs.get('threads', 96)
        
        conda_env = self.config.get_conda_env("metawrap")
        if not conda_env:
            raise ValueError("MetaWRAP conda environment not configured")
        
        # Build MetaWRAP binning command
        binner_flags = []
        for binner in binners:
            if binner.lower() in ['metabat1', 'metabat2', 'maxbin2', 'concoct']:
                binner_flags.append(f"--{binner.lower()}")
        
        commands = [
            f"metawrap binning -o {output_dir} -t {threads} -a {contigs_file} {' '.join(binner_flags)} {reads_file} --single-end --run-checkm"
        ]
        
        try:
            utils.execute_conda_commands(
                self.conda_script,
                conda_env,
                commands,
                timeout=kwargs.get('timeout')
            )
            
            logger.info(f"MetaWRAP binning completed: {output_dir}")
            return str(output_dir)
            
        except Exception as e:
            logger.error(f"MetaWRAP binning failed: {e}")
            raise

def run(args, config_manager: config.ConfigManager):
    """Main entry point for binning module"""
    
    if args.binning_tool == 'metabat2':
        binner = MetaBat2Binner(config_manager)
        bins_dir = binner.run(
            args.contigs,
            args.reads,
            args.output,
            threads=getattr(args, 'threads', None)
        )
    elif args.binning_tool == 'concoct':
        binner = ConcoTBinner(config_manager)
        bins_dir = binner.run(
            args.contigs,
            args.reads,
            args.output,
            threads=getattr(args, 'threads', None)
        )
    elif args.binning_tool == 'metawrap':
        binner = MetaWRAPBinner(config_manager)
        bins_dir = binner.run(
            args.contigs,
            args.reads,
            args.output,
            threads=getattr(args, 'threads', None)
        )
    else:
        raise ValueError(f"Unknown binning tool: {args.binning_tool}")
    
    logger.info(f"Binning completed successfully: {bins_dir}")
    return bins_dir