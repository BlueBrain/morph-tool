#!/usr/bin/env python

import imp
import sys
from pathlib import Path

from setuptools import find_packages, setup

# read the contents of the README file
readme_path = Path(__file__).resolve().parent / "README.rst"
with open(readme_path, encoding='utf-8') as f:
    long_description = f.read()

base_extras = [
    'bluepyopt>=1.9.37',
    'click>=6.7',
    'deprecation>=2.1.0',
    'more-itertools>=8.6.0',
    'morphio>=2.7.0',
    'neurom>=1.8.0',
    'numpy>=1.14',
    'pandas>=1.0.3',
    'xmltodict>=0.12.0',
]
plot_extras = ['plotly>=4.1.0']
parallel_extras = ['dask[bag]>=2.19.0']
nrn_extras = ['neuron>=7.8']

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
    install_requires=base_extras,
    extras_require={
        'nrn': nrn_extras,
        'parallel': parallel_extras,
        'plot': plot_extras,
        'all': plot_extras + parallel_extras + nrn_extras,
        'docs': ['sphinx', 'sphinx-bluebrain-theme'],
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
