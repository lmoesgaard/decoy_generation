"""
Setup script for decoy generation package.

This allows proper installation of the package with pip install -e .
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "A package for generating molecular decoys based on property similarity."

setup(
    name="decoy-generation",
    version="1.0.0",
    author="Decoy Generation Team",
    description="A package for generating molecular decoys based on property similarity",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "rdkit",
        "pandas",
        "numpy",
        "openbabel-wheel",  # Alternative to openbabel for easier installation
    ],
    entry_points={
        "console_scripts": [
            "decoy-generation=decoy_generation.main:main",
            "sample-props=decoy_generation.sample_props:main",
            "create-config=decoy_generation.write_input:create_example_configs",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)