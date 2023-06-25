from setuptools import setup, find_packages
import pathlib
import codecs
import os

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text(encoding='utf-8')

VERSION = '0.1.0'
DESCRIPTION = 'Processing open-source amplicon database for ecology'


setup(
    name='madt',
    version=VERSION,
    author="Yufei Zeng",
    author_email="yfzeng0827@hotmail.com",
    url='https://pypi.org/project/madt/',
    description=DESCRIPTION,
    long_description=README + '\n',
    long_description_content_type="text/markdown",
    license="GNU General Public License v3.0",
    classifiers=[ # see https://pypi.org/classifiers/
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords=['madt'],
    scripts=[
        'madt/madt_DataDownload.py', 
        'madt/madt_DataIndex.py', 
        'madt/madt_LibCheck.py', 
        'madt/madt_LibMap.py'
        ],
    # package_dir={'madt':'core'},
    # packages=find_packages(),
    packages=["madt"],
    include_package_data=True,
    install_requires=[
        'biopython',
        'pandas',
    ]
)
