#!/usr/bin/env python

import imp
import sys
from pathlib import Path

from setuptools import setup, find_packages

# read the contents of the README file
readme_path = Path(__file__).resolve().parent / "README.rst"
with open(readme_path, encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="morph-tool",
    author="Blue Brain Project, EPFL",
    description="A collection of CLIs and python function related to morphology handling",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/bluebrain/morph-tool",
    entry_points={
        'console_scripts': ['morph-tool=morph_tool.cli:cli']},
    license="LGPLv3",
    install_requires=[
        'click>=6.7',
        'functools32>=3.2;python_version<"3.0"',
        'morphio>=2.3.4',
        'numpy>=1.14',
        'neurom>=1.4.15',
        'pandas>=1.0.3',
    ],
    extras_require={
        'all': ['bluepyopt>=1.6',
                'neuron>=7.8',
                ],
    },
    python_requires='>=3.6',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
)
