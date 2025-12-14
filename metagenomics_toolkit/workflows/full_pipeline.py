import logging
from pathlib import Path
from typing import Optional, Dict, Any
import time

from ..core import utils, config
from ..modules import assembly, binning, classification, annotation, quality_control, utilities

logger = logging.getLogger(__name__)

class MetagenomicsPipeline:
    """Complete metagenomics analysis pipeline"""
    
    def __init__(self, config_manager: config.ConfigManager):
        self.config = config_manager
    
    def run_full_pipeline(self, input_file: str, output_dir: str, 
                         assembler: str = "flye", binner: str = "metabat2",
                         skip_qc: bool = False, **kwargs) -> Dict[str, str]:
        """Run complete metagenomics pipeline
        
        Args:
            input_file: Input FASTQ file
            output_dir: Base output directory
            assembler: Assembly tool to use
            binner: Binning tool to use
            skip_qc: Whether to skip quality control steps
            **kwargs: Additional parameters
            
        Returns:
            Dictionary with paths to all results
        """
        start_time = time.time()
        
        logger.info("=" * 60)
        logger.info("STARTING METAGENOMICS PIPELINE")
        logger.info("=" * 60)
        logger.info(f"Input: {input_file}")
        logger.info(f"Output: {output_dir}")
        logger.info(f"Assembler: {assembler}")
        logger.info(f"Binner: {binner}")
        logger.info(f"Skip QC: {skip_qc}")
        
        # Validate inputs
        utils.validate_file_exists(input_file, "Input FASTQ file")
        base_output = Path(output_dir)
        utils.ensure_directory(base_output)
        
        results = {
            'input_file': input_file,
            'output_directory': str(base_output)
        }
        
        try:
            # Step 1: Assembly
            logger.info("\n" + "="*40)
            logger.info("STEP 1: ASSEMBLY")
            logger.info("="*40)
            
            assembly_dir = base_output / "01_assembly"
            assembly_file = self._run_assembly(
                input_file, str(assembly_dir), assembler, **kwargs
            )
            results['assembly'] = assembly_file
            
            # Step 1.5: Assembly statistics
            logger.info("\n" + "-"*30)
            logger.info("Calculating assembly statistics")
            logger.info("-"*30)
            
            stats_calc = utilities.AssemblyStatsCalculator()
            assembly_stats = stats_calc.calculate_stats(assembly_file)
            stats_file = base_output / "01_assembly" / "assembly_stats.txt"
            
            with open(stats_file, 'w') as f:
                f.write(stats_calc.format_stats(assembly_stats))
            
            logger.info(f"Assembly statistics saved: {stats_file}")
            results['assembly_stats'] = str(stats_file)
            
            # Step 2: Binning
            logger.info("\n" + "="*40)
            logger.info("STEP 2: BINNING")
            logger.info("="*40)
            
            binning_dir = base_output / "02_binning"
            bins_dir = self._run_binning(
                assembly_file, input_file, str(binning_dir), binner, **kwargs
            )
            results['binning'] = bins_dir
            
            # Step 3: Quality Control (optional)
            if not skip_qc:
                logger.info("\n" + "="*40)
                logger.info("STEP 3: QUALITY CONTROL")
                logger.info("="*40)
                
                qc_dir = base_output / "03_quality_control"
                qc_results = self._run_quality_control(
                    bins_dir, str(qc_dir), **kwargs
                )
                results['quality_control'] = qc_results
            
            # Step 4: Classification
            logger.info("\n" + "="*40)
            logger.info("STEP 4: TAXONOMIC CLASSIFICATION")
            logger.info("="*40)
            
            classification_dir = base_output / "04_classification"
            classification_results = self._run_classification(
                bins_dir, str(classification_dir), **kwargs
            )
            results['classification'] = classification_results
            
            # Step 5: Annotation
            logger.info("\n" + "="*40)
            logger.info("STEP 5: GENOME ANNOTATION")
            logger.info("="*40)
            
            annotation_dir = base_output / "05_annotation"
            annotation_results = self._run_annotation(
                bins_dir, str(annotation_dir), **kwargs
            )
            results['annotation'] = annotation_results
            
            # Generate final summary
            self._generate_summary(results, base_output)
            
            end_time = time.time()
            runtime = end_time - start_time
            
            logger.info("\n" + "="*60)
            logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            logger.info("="*60)
            logger.info(f"Total runtime: {runtime:.2f} seconds ({runtime/3600:.2f} hours)")
            logger.info(f"Results directory: {base_output}")
            
            return results
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
    
    def _run_assembly(self, input_file: str, output_dir: str, 
                     assembler: str, **kwargs) -> str:
        """Run assembly step"""
        
        # Create mock args object for assembly module
        class AssemblyArgs:
            def __init__(self):
                self.assembly_tool = assembler
                self.input = input_file
                self.output = output_dir
                self.threads = kwargs.get('threads', 4)
                if assembler == 'flye':
                    self.read_type = kwargs.get('read_type', 'nano-raw')
        
        args = AssemblyArgs()
        return assembly.run(args, self.config)
    
    def _run_binning(self, contigs_file: str, reads_file: str, 
                    output_dir: str, binner: str, **kwargs) -> str:
        """Run binning step"""
        
        # Create mock args object for binning module
        class BinningArgs:
            def __init__(self):
                self.binning_tool = binner
                self.contigs = contigs_file
                self.reads = reads_file
                self.output = output_dir
                self.threads = kwargs.get('threads', 8)
        
        args = BinningArgs()
        return binning.run(args, self.config)
    
    def _run_quality_control(self, bins_dir: str, output_dir: str, 
                           **kwargs) -> str:
        """Run quality control step"""
        
        # Create mock args object for QC module
        class QCArgs:
            def __init__(self):
                self.qc_tool = 'checkm'
                self.bins_dir = bins_dir
                self.output = output_dir
                self.extension = kwargs.get('bin_extension', 'fa')
                self.threads = kwargs.get('threads', 8)
        
        args = QCArgs()
        return quality_control.run(args, self.config)
    
    def _run_classification(self, bins_dir: str, output_dir: str, 
                          **kwargs) -> str:
        """Run classification step"""
        
        # Create mock args object for classification module
        class ClassificationArgs:
            def __init__(self):
                self.classify_tool = 'gtdbtk'
                self.bins_dir = bins_dir
                self.output = output_dir
                self.extension = kwargs.get('bin_extension', 'fa')
                self.threads = kwargs.get('threads', 8)
        
        args = ClassificationArgs()
        return classification.run(args, self.config)
    
    def _run_annotation(self, bins_dir: str, output_dir: str, 
                       **kwargs) -> str:
        """Run annotation step"""
        
        # Create mock args object for annotation module
        class AnnotationArgs:
            def __init__(self):
                self.annotate_tool = 'prokka'
                self.input = bins_dir
                self.output = output_dir
                self.continue_on_error = True
        
        args = AnnotationArgs()
        return annotation.run(args, self.config)
    
    def _generate_summary(self, results: Dict[str, str], output_dir: Path):
        """Generate pipeline summary report"""
        
        summary_file = output_dir / "PIPELINE_SUMMARY.txt"
        
        with open(summary_file, 'w') as f:
            f.write("METAGENOMICS PIPELINE SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Input file: {results['input_file']}\n")
            f.write(f"Output directory: {results['output_directory']}\n\n")
            
            f.write("PIPELINE STEPS COMPLETED:\n")
            f.write("-" * 30 + "\n")
            
            if 'assembly' in results:
                f.write(f"✓ Assembly: {results['assembly']}\n")
            
            if 'assembly_stats' in results:
                f.write(f"✓ Assembly Statistics: {results['assembly_stats']}\n")
            
            if 'binning' in results:
                f.write(f"✓ Binning: {results['binning']}\n")
            
            if 'quality_control' in results:
                f.write(f"✓ Quality Control: {results['quality_control']}\n")
            
            if 'classification' in results:
                f.write(f"✓ Taxonomic Classification: {results['classification']}\n")
            
            if 'annotation' in results:
                f.write(f"✓ Genome Annotation: {results['annotation']}\n")
            
            f.write("\nPIPELINE COMPLETED SUCCESSFULLY!\n")
        
        logger.info(f"Pipeline summary saved: {summary_file}")

def run(args, config_manager: config.ConfigManager):
    """Entry point for full pipeline workflow"""
    
    pipeline = MetagenomicsPipeline(config_manager)
    
    results = pipeline.run_full_pipeline(
        input_file=args.input,
        output_dir=args.output,
        assembler=args.assembler,
        binner=args.binning,
        skip_qc=args.skip_qc,
        threads=getattr(args, 'threads', 4),
        read_type=getattr(args, 'read_type', 'nano-raw'),
        bin_extension=getattr(args, 'bin_extension', 'fa')
    )
    
    return results