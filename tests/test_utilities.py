import pytest
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metagenomics_toolkit.modules.utilities import (
    AssemblyStatsCalculator,
    SequenceRenamer,
    GeneFinder
)

class TestAssemblyStatsCalculator:
    """Test assembly statistics calculation"""
    
    def create_test_fasta(self, sequences):
        """Create a test FASTA file with given sequences"""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        
        records = []
        for i, seq_str in enumerate(sequences):
            record = SeqRecord(Seq(seq_str), id=f"contig_{i+1}", description=f"test contig {i+1}")
            records.append(record)
        
        SeqIO.write(records, temp_file, "fasta")
        temp_file.close()
        return temp_file.name
    
    def test_basic_stats(self):
        """Test basic assembly statistics"""
        # Create test sequences
        sequences = [
            "ATCGATCGATCGATCG",  # 16 bp
            "GCTAGCTAGCTA",       # 12 bp
            "TTTTAAAA",           # 8 bp
            "AAAAAAAAAAAAAAAA"    # 16 bp
        ]
        
        fasta_file = self.create_test_fasta(sequences)
        
        try:
            calculator = AssemblyStatsCalculator()
            stats = calculator.calculate_stats(fasta_file)
            
            # Verify basic stats
            assert stats['num_contigs'] == 4
            assert stats['total_length'] == 52  # 16 + 12 + 8 + 16
            assert stats['min_length'] == 8
            assert stats['max_length'] == 16
            assert stats['mean_length'] == 13.0  # 52 / 4
            
            # N50 should be 16 (largest sequence in top 50%)
            assert stats['n50'] == 16
            
        finally:
            Path(fasta_file).unlink()
    
    def test_gc_content(self):
        """Test GC content calculation"""
        sequences = [
            "GGGGCCCC",  # 100% GC
            "AAAATTTT"   # 0% GC
        ]
        
        fasta_file = self.create_test_fasta(sequences)
        
        try:
            calculator = AssemblyStatsCalculator()
            stats = calculator.calculate_stats(fasta_file)
            
            # Should be 50% GC (8 GC bases out of 16 total)
            assert abs(stats['gc_content'] - 50.0) < 0.1
            
        finally:
            Path(fasta_file).unlink()
    
    def test_format_stats(self):
        """Test statistics formatting"""
        stats = {
            'file': 'test.fasta',
            'num_contigs': 100,
            'total_length': 1000000,
            'min_length': 500,
            'max_length': 50000,
            'mean_length': 10000.0,
            'n50': 15000,
            'n90': 2000,
            'gc_content': 45.5
        }
        
        calculator = AssemblyStatsCalculator()
        formatted = calculator.format_stats(stats)
        
        assert "test.fasta" in formatted
        assert "100" in formatted
        assert "1,000,000" in formatted  # Check formatting with commas
        assert "45.50%" in formatted  # Two decimal places expected

class TestSequenceRenamer:
    """Test sequence renaming utilities"""
    
    def create_test_fasta(self, sequences_with_names):
        """Create test FASTA file with specific sequence names"""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        
        records = []
        for seq_id, seq_str in sequences_with_names:
            record = SeqRecord(Seq(seq_str), id=seq_id, description=f"description for {seq_id}")
            records.append(record)
        
        SeqIO.write(records, temp_file, "fasta")
        temp_file.close()
        return temp_file.name
    
    def test_rename_single_file(self):
        """Test renaming sequences in a single file"""
        sequences = [
            ("contig_1", "ATCGATCG"),
            ("contig_2", "GCTAGCTA")
        ]
        
        fasta_file = self.create_test_fasta(sequences)
        
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                renamer = SequenceRenamer()
                output_dir = renamer.rename_sequences(fasta_file, "sample1", tmpdir)
                
                # Check output file exists
                output_files = list(Path(output_dir).glob("*.fasta"))
                assert len(output_files) == 1
                
                # Check sequences were renamed
                renamed_records = list(SeqIO.parse(output_files[0], "fasta"))
                assert len(renamed_records) == 2
                assert renamed_records[0].id.startswith("sample1_")
                assert renamed_records[1].id.startswith("sample1_")
                
        finally:
            Path(fasta_file).unlink()

class TestGeneFinder:
    """Test gene finding utilities"""
    
    def create_test_fasta_with_genes(self):
        """Create test FASTA with gene annotations"""
        sequences = [
            ("gene1_rpoB_function", "ATCGATCGATCG"),
            ("gene2_16S_ribosomal", "GCTAGCTAGCTA"),
            ("gene3_unknown", "TTTTAAAA"),
            ("gene4_rpoB_variant", "GGGGCCCC")
        ]
        
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        
        records = []
        for seq_id, seq_str in sequences:
            record = SeqRecord(Seq(seq_str), id=seq_id, description=f"description for {seq_id}")
            records.append(record)
        
        SeqIO.write(records, temp_file, "fasta")
        temp_file.close()
        return temp_file.name
    
    def test_find_genes_case_insensitive(self):
        """Test gene finding (case insensitive)"""
        fasta_file = self.create_test_fasta_with_genes()
        
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as output_file:
                finder = GeneFinder()
                num_found = finder.find_genes(
                    fasta_file, 
                    ["rpoB", "16S"], 
                    output_file.name,
                    case_sensitive=False
                )
                
                assert num_found == 3  # 2 rpoB + 1 16S
                
                # Check output content
                with open(output_file.name) as f:
                    content = f.read()
                    assert "rpoB" in content
                    assert "16S" in content
                    assert "gene1_rpoB_function" in content
                    assert "gene4_rpoB_variant" in content
                    
        finally:
            Path(fasta_file).unlink()
            Path(output_file.name).unlink()
    
    def test_find_genes_case_sensitive(self):
        """Test gene finding (case sensitive)"""
        fasta_file = self.create_test_fasta_with_genes()
        
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as output_file:
                finder = GeneFinder()
                num_found = finder.find_genes(
                    fasta_file, 
                    ["rpoB"],  # Exact case match
                    output_file.name,
                    case_sensitive=True
                )
                
                assert num_found == 2  # 2 exact rpoB matches
                
        finally:
            Path(fasta_file).unlink()
            Path(output_file.name).unlink()