# Metagenomics Toolkit

A comprehensive, modular Python wrapper for metagenomic analysis pipelines, integrating assembly, binning, classification, annotation, and quality control tools.


## Scientific Publications Applying this Pipeline
Metagenomes and Metagenome-Assembled Genomes from Microbial Communities in a Biological Nutrient Removal Plant Operated at Los Angeles County Sanitation District (LACSD) with High and Low Dissolved Oxygen Conditions. doi: https://doi.org/10.1101/2025.11.10.687646. 

Metagenomes and Metagenome-Assembled Genomes from Microbial Communities in a Biological Nutrient Removal Plant Operated at Hamptons Road Sanitation District (HRSD) with High and Low Dissolved Oxygen Conditions. doi: https://doi.org/10.1101/2025.11.10.687637

## Features

- **Modular Design**: Each analysis step is a separate module that can be run independently
- **Multiple Tool Support**: Supports various tools for each step (FLYE, SPAdes, MEGAHIT for assembly; MetaBAT2, CONCOCT for binning; etc.)
- **Command-Line Interface**: Easy-to-use CLI with comprehensive help and options
- **Configuration Management**: Flexible configuration system for tool paths and parameters
- **Error Handling**: Robust error handling with informative error messages
- **Logging**: Comprehensive logging for troubleshooting and monitoring
- **Conda Integration**: Seamless integration with conda environments

## Supported Tools

### Assembly
- **FLYE**: Long-read metagenomic assembly
- **SPAdes**: Short-read metagenomic assembly  
- **MEGAHIT**: Fast and memory-efficient metagenome assembly

### Binning
- **MetaBAT2**: Efficient binning of large metagenomic datasets
- **CONCOCT**: Clustering contigs on coverage and composition
- **MetaWRAP**: Comprehensive binning wrapper

### Classification
- **GTDB-tk**: Taxonomic classification using GTDB

### Annotation
- **Prokka**: Genome annotation for prokaryotes

### Quality Control
- **CheckM**: Assess genome quality and contamination

### Utilities
- Assembly statistics calculation
- Sequence renaming with prefixes
- Gene sequence finding and extraction

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/your-username/metagenomics-toolkit.git
cd metagenomics-toolkit

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Prerequisites

- Python 3.8+
- Conda/Miniconda for bioinformatics tools
- Bioinformatics tools (FLYE, MetaBAT2, GTDB-tk, etc.) installed in conda environments

## Quick Start

### Full Pipeline
Run a complete metagenomics pipeline:

```bash
metagenomics-toolkit pipeline \
    --input reads.fastq \
    --assembler flye \
    --binning metabat2 \
    --output results/ \
    --threads 8
```

### Individual Steps

#### Assembly
```bash
# FLYE assembly
metagenomics-toolkit assembly flye \
    --input reads.fastq \
    --output assembly/ \
    --read-type nano-raw

# SPAdes assembly  
metagenomics-toolkit assembly spades \
    --input read1.fastq,read2.fastq \
    --output assembly/
```

#### Binning
```bash
# MetaBAT2 binning
metagenomics-toolkit binning metabat2 \
    --contigs contigs.fasta \
    --reads reads.fastq \
    --output bins/
```

#### Classification
```bash
# GTDB-tk classification
metagenomics-toolkit classify gtdbtk \
    --bins-dir bins/ \
    --output classification/
```

#### Annotation
```bash
# Prokka annotation
metagenomics-toolkit annotate prokka \
    --input bins/ \
    --output annotations/
```

#### Quality Control
```bash
# CheckM assessment
metagenomics-toolkit qc checkm \
    --bins-dir bins/ \
    --output qc_results/
```

#### Utilities
```bash
# Calculate assembly statistics
metagenomics-toolkit utils stats \
    --input assembly.fasta \
    --output stats.txt

# Rename sequences with prefix
metagenomics-toolkit utils rename \
    --input sequences.fasta \
    --prefix sample1 \
    --output renamed/

# Find gene sequences
metagenomics-toolkit utils find-genes \
    --fasta annotations.fasta \
    --genes rpoB 16S recA \
    --output found_genes.txt
```

## Configuration

Create a custom configuration file to specify tool paths and parameters:

```json
{
  "conda": {
    "activate_script": "/opt/conda/etc/profile.d/conda.sh",
    "environments": {
      "flye": "/opt/conda/envs/flye",
      "metabat2": "/opt/conda/envs/metabat2",
      "gtdbtk": "/opt/conda/envs/gtdbtk"
    }
  },
  "tools": {
    "flye": {
      "default_threads": 20,
      "meta_mode": true
    },
    "metabat2": {
      "min_contig_length": 2000
    }
  }
}
```

Use with:
```bash
metagenomics-toolkit --config config.json pipeline --input reads.fastq --output results/
```

## Output Structure

The pipeline creates a structured output directory:

```
results/
├── 01_assembly/
│   ├── assembly.fasta
│   └── assembly_stats.txt
├── 02_binning/
│   └── bins/
├── 03_quality_control/
│   └── checkm_results/
├── 04_classification/
│   └── gtdbtk_results/
├── 05_annotation/
│   └── prokka_results/
└── PIPELINE_SUMMARY.txt
```

## Requirements

### Python Dependencies
- biopython>=1.79
- click>=8.0.0
- pyyaml>=6.0
- psutil>=5.8.0

### Bioinformatics Tools
All bioinformatics tools should be installed in separate conda environments:

- FLYE (for assembly)
- SPAdes (for assembly)
- MEGAHIT (for assembly)
- MetaBAT2 (for binning)
- MetaWRAP (for binning)
- GTDB-tk (for classification)
- Prokka (for annotation)
- CheckM (for quality control)
- minimap2, samtools (for read mapping)

## Development

### Running Tests

```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Run tests with coverage
pytest tests/ --cov=metagenomics_toolkit
```

### Code Formatting

```bash
# Format code with black
black metagenomics_toolkit/

# Lint code with flake8
flake8 metagenomics_toolkit/
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```
Enuh, B.M. (2024). Metagenomics Toolkit: A comprehensive Python wrapper for metagenomic analysis pipelines. 
```

## Support

For questions, bug reports, or feature requests, please open an issue on GitHub.

## Acknowledgments

- Built on top of excellent bioinformatics tools from the community
- Inspired by the need for reproducible and standardized metagenomic workflows
- Thanks to the developers of FLYE, MetaBAT2, GTDB-tk, Prokka, CheckM, and other tools
