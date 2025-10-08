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
    author="Lars Moesgaard",
    author_email="your.email@example.com",  # Update with your email
    description="Molecular decoy generation with property-based filtering and structural similarity analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lmoesgaard/decoy_generation",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "rdkit>=2023.3.1",
        "pandas>=1.5.0",
        "numpy>=1.21.0",
        "matplotlib>=3.5.0",
        "openbabel-wheel>=3.1.0",
        "scipy>=1.9.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "black>=22.0.0",
            "flake8>=5.0.0",
            "mypy>=0.900",
        ],
        "docs": [
            "sphinx>=5.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "decoy-generation=main:main",
            "sample-props=sample_props:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="chemistry, molecular, decoys, drug-discovery, cheminformatics",
    project_urls={
        "Bug Reports": "https://github.com/lmoesgaard/decoy_generation/issues",
        "Source": "https://github.com/lmoesgaard/decoy_generation",
        "Documentation": "https://github.com/lmoesgaard/decoy_generation#readme",
    },
    include_package_data=True,
    package_data={
        "decoy_generation": [
            "examples/*.json",
            "examples/*.md",
        ],
    },
)

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
        "matplotlib",
        "scipy",
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