#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

VERSION = imp.load_source("", "morph_tool/version.py").__version__

PLOTLY_EXTRAS = [
    'plotly>=4.1.0',
    'pandas>=0.24',
    'bluepy>=0.14',
]

setup(
    name="morph-tool",
    author="BlueBrain NSE",
    author_email="bbp-ou-nse@groupes.epfl.ch",
    version=VERSION,
    description="A collection of CLIs and python function related to morphology handling",
    url="https://github.com/bluebrain/morph-tool",
    entry_points={
        'console_scripts': ['morph-tool=morph_tool.cli:cli']},
    license="BBP-internal-confidential",
    install_requires=[
        'click>=6.7',
        'functools32>=3.2;python_version<"3.0"',
        'morphio>=2.3.4',
        'numpy>=1.14',
        'neurom>=1.4.15',
    ],
    extras_require={
        'all': ['bluepyopt>=1.6'] + PLOTLY_EXTRAS,
        'plotly': PLOTLY_EXTRAS,
    },
    python_requires='>=3.6',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)
