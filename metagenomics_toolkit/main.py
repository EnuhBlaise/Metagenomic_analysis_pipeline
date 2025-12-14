#!/usr/bin/env python3

import argparse
import sys
import os
import logging
from pathlib import Path

from .core import config, logger
from .modules import assembly, binning, classification, annotation, quality_control, utilities

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logger.setup_logger(level)

def create_parser():
    parser = argparse.ArgumentParser(
        prog='metagenomics-toolkit',
        description='Comprehensive Metagenomics Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete pipeline with FLYE assembly
  metagenomics-toolkit pipeline --input reads.fastq --assembler flye --output results/

  # Assembly only
  metagenomics-toolkit assembly flye --input reads.fastq --output assembly/

  # Binning with MetaBAT2  
  metagenomics-toolkit binning metabat2 --contigs contigs.fasta --reads reads.fastq --output bins/

  # Classification with GTDB-tk
  metagenomics-toolkit classify gtdbtk --bins-dir bins/ --output classification/
        """
    )
    
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('--output-dir', help='Base output directory')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Pipeline command
    pipeline_parser = subparsers.add_parser('pipeline', help='Run complete metagenomics pipeline')
    pipeline_parser.add_argument('--input', required=True, help='Input FASTQ file or directory')
    pipeline_parser.add_argument('--assembler', choices=['flye', 'spades', 'megahit', 'metamdbg'], 
                                 default='flye', help='Assembly tool')
    pipeline_parser.add_argument('--binning', choices=['metabat2', 'concoct', 'metawrap'], 
                                 default='metabat2', help='Binning tool')
    pipeline_parser.add_argument('--output', required=True, help='Output directory')
    pipeline_parser.add_argument('--skip-qc', action='store_true', help='Skip quality control steps')
    
    # Assembly subcommands
    assembly_parser = subparsers.add_parser('assembly', help='Assembly operations')
    assembly_subparsers = assembly_parser.add_subparsers(dest='assembly_tool')
    
    for tool in ['flye', 'spades', 'megahit', 'metamdbg']:
        tool_parser = assembly_subparsers.add_parser(tool, help=f'{tool.upper()} assembly')
        tool_parser.add_argument('--input', required=True, help='Input FASTQ file')
        tool_parser.add_argument('--output', required=True, help='Output directory')
        if tool == 'flye':
            tool_parser.add_argument('--read-type', choices=['nano-raw', 'nano-corr', 'pacbio'], 
                                   default='nano-raw', help='Read type for FLYE')
    
    # Binning subcommands
    binning_parser = subparsers.add_parser('binning', help='Binning operations')
    binning_subparsers = binning_parser.add_subparsers(dest='binning_tool')
    
    for tool in ['metabat2', 'concoct', 'metawrap']:
        tool_parser = binning_subparsers.add_parser(tool, help=f'{tool.upper()} binning')
        tool_parser.add_argument('--contigs', required=True, help='Input contigs file')
        tool_parser.add_argument('--reads', required=True, help='Input reads file')
        tool_parser.add_argument('--output', required=True, help='Output directory')
    
    # Classification subcommands
    classify_parser = subparsers.add_parser('classify', help='Taxonomic classification')
    classify_subparsers = classify_parser.add_subparsers(dest='classify_tool')
    
    gtdbtk_parser = classify_subparsers.add_parser('gtdbtk', help='GTDB-tk classification')
    gtdbtk_parser.add_argument('--bins-dir', required=True, help='Directory containing bins')
    gtdbtk_parser.add_argument('--output', required=True, help='Output directory')
    gtdbtk_parser.add_argument('--extension', default='fa', help='File extension for bins')
    
    # Annotation subcommands
    annotate_parser = subparsers.add_parser('annotate', help='Genome annotation')
    annotate_subparsers = annotate_parser.add_subparsers(dest='annotate_tool')
    
    prokka_parser = annotate_subparsers.add_parser('prokka', help='Prokka annotation')
    prokka_parser.add_argument('--input', required=True, help='Input genome file or directory')
    prokka_parser.add_argument('--output', required=True, help='Output directory')
    
    # Quality control subcommands
    qc_parser = subparsers.add_parser('qc', help='Quality control')
    qc_subparsers = qc_parser.add_subparsers(dest='qc_tool')
    
    checkm_parser = qc_subparsers.add_parser('checkm', help='CheckM quality assessment')
    checkm_parser.add_argument('--bins-dir', required=True, help='Directory containing bins')
    checkm_parser.add_argument('--output', required=True, help='Output directory')
    checkm_parser.add_argument('--extension', default='fa', help='File extension for bins')
    
    # Utilities subcommands
    utils_parser = subparsers.add_parser('utils', help='Utility operations')
    utils_subparsers = utils_parser.add_subparsers(dest='utils_tool')
    
    stats_parser = utils_subparsers.add_parser('stats', help='Calculate assembly statistics')
    stats_parser.add_argument('--input', required=True, help='Input FASTA file')
    stats_parser.add_argument('--output', help='Output file (default: stdout)')
    
    rename_parser = utils_subparsers.add_parser('rename', help='Rename sequences with prefix')
    rename_parser.add_argument('--input', required=True, help='Input file/directory')
    rename_parser.add_argument('--prefix', required=True, help='Prefix to add')
    rename_parser.add_argument('--output', help='Output directory')
    
    find_genes_parser = utils_subparsers.add_parser('find-genes', help='Find gene sequences')
    find_genes_parser.add_argument('--fasta', required=True, help='Input FASTA file')
    find_genes_parser.add_argument('--genes', nargs='+', required=True, help='Gene names to search')
    find_genes_parser.add_argument('--output', required=True, help='Output file')
    
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    setup_logging(args.verbose)
    
    # Load configuration
    config_manager = config.ConfigManager(args.config)
    
    try:
        # Route to appropriate module
        if args.command == 'pipeline':
            from .workflows import full_pipeline
            full_pipeline.run(args, config_manager)
            
        elif args.command == 'assembly':
            assembly.run(args, config_manager)
            
        elif args.command == 'binning':
            binning.run(args, config_manager)
            
        elif args.command == 'classify':
            classification.run(args, config_manager)
            
        elif args.command == 'annotate':
            annotation.run(args, config_manager)
            
        elif args.command == 'qc':
            quality_control.run(args, config_manager)
            
        elif args.command == 'utils':
            utilities.run(args, config_manager)
            
        else:
            parser.print_help()
            sys.exit(1)
            
    except KeyboardInterrupt:
        logging.error("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()