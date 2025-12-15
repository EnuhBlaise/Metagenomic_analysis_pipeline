#!/usr/bin/env python3

from setuptools import setup, find_packages
from pathlib import Path

# Read README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text() if (this_directory / "README.md").exists() else ""

# Read requirements
def read_requirements(filename):
    """Read requirements from file"""
    requirements_path = this_directory / filename
    if requirements_path.exists():
        with open(requirements_path) as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="metagenomics-toolkit",
    version="1.0.0",
    author="Blaise Manga Enuh",
    author_email="blaisonyl@gmail.com",
    description="A comprehensive modular metagenomics analysis pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EnuhBlaise/metagenomics-toolkit",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements("requirements.txt"),
    extras_require={
        "dev": read_requirements("requirements-dev.txt"),
        "test": read_requirements("requirements-test.txt"),
    },
    entry_points={
        "console_scripts": [
            "metagenomics-toolkit=metagenomics_toolkit.main:main",
            "meta-toolkit=metagenomics_toolkit.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "metagenomics_toolkit": [
            "configs/*.json",
            "scripts/*.sh",
        ],
    },
    zip_safe=False,
    keywords="metagenomics bioinformatics assembly binning classification annotation",
    project_urls={
        "Bug Reports": "https://github.com/EnuhBlaise/metagenomic_analysis_pipeline/issues",
        "Source": "https://github.com/EnuhBlaise/metagenomic_analysis_pipeline",
        "Documentation": "https://github.com/EnuhBlaise/metagenomic_analysis_pipeline/wiki",
    },
)
