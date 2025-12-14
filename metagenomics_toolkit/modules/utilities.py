import logging
from pathlib import Path
from typing import List, Dict, Any
import shutil
from Bio import SeqIO
import argparse

from ..core import utils, config

logger = logging.getLogger(__name__)

class AssemblyStatsCalculator:
    """Calculate assembly statistics"""
    
    @staticmethod
    def calculate_stats(fasta_file: str) -> Dict[str, Any]:
        """Calculate assembly statistics from FASTA file
        
        Args:
            fasta_file: Input FASTA file
            
        Returns:
            Dictionary with assembly statistics
        """
        utils.validate_file_exists(fasta_file, "FASTA file")
        
        logger.info(f"Calculating assembly statistics: {fasta_file}")
        
        sequences = []
        total_length = 0
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_len = len(record.seq)
                sequences.append(seq_len)
                total_length += seq_len
        except Exception as e:
            raise ValueError(f"Error reading FASTA file {fasta_file}: {e}")
        
        if not sequences:
            raise ValueError(f"No sequences found in {fasta_file}")
        
        # Sort sequences by length (descending)
        sequences.sort(reverse=True)
        
        num_contigs = len(sequences)
        min_length = min(sequences)
        max_length = max(sequences)
        mean_length = total_length / num_contigs
        
        # Calculate N50
        target_length = total_length * 0.5
        cumulative_length = 0
        n50 = 0
        
        for length in sequences:
            cumulative_length += length
            if cumulative_length >= target_length:
                n50 = length
                break
        
        # Calculate N90
        target_length_90 = total_length * 0.9
        cumulative_length = 0
        n90 = 0
        
        for length in sequences:
            cumulative_length += length
            if cumulative_length >= target_length_90:
                n90 = length
                break
        
        # Calculate GC content
        total_gc = 0
        total_bases = 0
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_str = str(record.seq).upper()
                gc_count = seq_str.count('G') + seq_str.count('C')
                total_gc += gc_count
                total_bases += len(seq_str)
        except Exception as e:
            logger.warning(f"Could not calculate GC content: {e}")
            gc_content = None
        else:
            gc_content = (total_gc / total_bases) * 100 if total_bases > 0 else 0
        
        stats = {
            'file': fasta_file,
            'num_contigs': num_contigs,
            'total_length': total_length,
            'min_length': min_length,
            'max_length': max_length,
            'mean_length': mean_length,
            'n50': n50,
            'n90': n90,
            'gc_content': gc_content
        }
        
        return stats
    
    @staticmethod
    def format_stats(stats: Dict[str, Any]) -> str:
        """Format statistics for display
        
        Args:
            stats: Statistics dictionary
            
        Returns:
            Formatted statistics string
        """
        lines = [
            f"Assembly Statistics for: {Path(stats['file']).name}",
            "=" * 50,
            f"Number of contigs: {stats['num_contigs']:,}",
            f"Total length: {stats['total_length']:,} bp",
            f"Min contig length: {stats['min_length']:,} bp",
            f"Max contig length: {stats['max_length']:,} bp",
            f"Mean contig length: {stats['mean_length']:.2f} bp",
            f"N50: {stats['n50']:,} bp",
            f"N90: {stats['n90']:,} bp"
        ]
        
        if stats['gc_content'] is not None:
            lines.append(f"GC content: {stats['gc_content']:.2f}%")
        
        return "\n".join(lines)

class SequenceRenamer:
    """Rename sequences with prefixes"""
    
    @staticmethod
    def rename_sequences(input_path: str, prefix: str, output_dir: str = None) -> str:
        """Rename sequences in FASTA files with prefix
        
        Args:
            input_path: Input file or directory
            prefix: Prefix to add to sequence names
            output_dir: Output directory (optional)
            
        Returns:
            Path to output directory or file
        """
        input_path_obj = Path(input_path)
        
        if input_path_obj.is_file():
            files_to_process = [input_path_obj]
        elif input_path_obj.is_dir():
            # Find FASTA files
            extensions = ['fa', 'fasta', 'fna']
            files_to_process = []
            for ext in extensions:
                files_to_process.extend(utils.find_files_with_extension(input_path, ext))
        else:
            raise ValueError(f"Input path does not exist: {input_path}")
        
        if not files_to_process:
            raise ValueError("No FASTA files found to process")
        
        # Determine output directory
        if output_dir is None:
            if input_path_obj.is_file():
                output_dir = input_path_obj.parent / f"{input_path_obj.stem}_renamed"
            else:
                output_dir = input_path_obj.parent / f"{input_path_obj.name}_renamed"
        
        output_path = Path(output_dir)
        utils.ensure_directory(output_path)
        
        logger.info(f"Renaming {len(files_to_process)} files with prefix '{prefix}'")
        
        for file_path in files_to_process:
            output_file = output_path / f"{prefix}_{file_path.name}"
            
            try:
                with open(output_file, 'w') as out_handle:
                    for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
                        new_id = f"{prefix}_{record.id}"
                        record.id = new_id
                        record.description = f"{prefix}_{record.description}"
                        SeqIO.write(record, out_handle, "fasta")
                
                logger.info(f"Renamed sequences in {file_path.name} -> {output_file.name}")
                
            except Exception as e:
                logger.error(f"Error processing {file_path}: {e}")
                raise
        
        return str(output_path)

class GeneFinder:
    """Find gene sequences in FASTA files"""
    
    @staticmethod
    def find_genes(fasta_file: str, gene_names: List[str], 
                   output_file: str, case_sensitive: bool = False) -> int:
        """Find and extract gene sequences
        
        Args:
            fasta_file: Input FASTA file
            gene_names: List of gene names to search for
            output_file: Output file path
            case_sensitive: Whether search should be case sensitive
            
        Returns:
            Number of sequences found
        """
        utils.validate_file_exists(fasta_file, "FASTA file")
        
        logger.info(f"Searching for {len(gene_names)} genes in {fasta_file}")
        
        sequences_found = []
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = record.id
                sequence = str(record.seq)
                
                # Check if any gene name is in the record header
                for gene_name in gene_names:
                    if case_sensitive:
                        if gene_name in header:
                            match = True
                    else:
                        if gene_name.lower() in header.lower():
                            match = True
                    
                    if 'match' in locals() and match:
                        sequences_found.append({
                            'gene_name': gene_name,
                            'header': header,
                            'sequence': sequence,
                        })
                        del match
                        break
        
        except Exception as e:
            raise ValueError(f"Error reading FASTA file {fasta_file}: {e}")
        
        # Write results
        utils.ensure_directory(Path(output_file).parent)
        
        with open(output_file, 'w') as f:
            for info in sequences_found:
                f.write(f"Gene Name: {info['gene_name']}\n")
                f.write(f"Header: {info['header']}\n")
                f.write(f"Sequence: {info['sequence']}\n\n")
        
        logger.info(f"Found {len(sequences_found)} gene sequences, saved to {output_file}")
        return len(sequences_found)

def run(args, config_manager: config.ConfigManager):
    """Main entry point for utilities module"""
    
    if args.utils_tool == 'stats':
        # Calculate assembly statistics
        calculator = AssemblyStatsCalculator()
        stats = calculator.calculate_stats(args.input)
        formatted_stats = calculator.format_stats(stats)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(formatted_stats)
            logger.info(f"Statistics saved to {args.output}")
        else:
            print(formatted_stats)
        
        return args.output or "stdout"
    
    elif args.utils_tool == 'rename':
        # Rename sequences
        renamer = SequenceRenamer()
        output_dir = renamer.rename_sequences(
            args.input,
            args.prefix,
            args.output
        )
        logger.info(f"Sequences renamed successfully: {output_dir}")
        return output_dir
    
    elif args.utils_tool == 'find-genes':
        # Find gene sequences
        finder = GeneFinder()
        num_found = finder.find_genes(
            args.fasta,
            args.genes,
            args.output,
            case_sensitive=getattr(args, 'case_sensitive', False)
        )
        logger.info(f"Found {num_found} gene sequences")
        return args.output
    
    else:
        raise ValueError(f"Unknown utility tool: {args.utils_tool}")