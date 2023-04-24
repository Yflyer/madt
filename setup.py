from setuptools import setup, find_packages

setup(
    name='madt',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'madt_DataIndex = madt.madt_DataIndex:main',
            'madt_LibCheck = madt.madt_LibCheck:main',
            'madt_LipMap = madt.madt_LipMap:main'
        ]
    },
    install_requires=[
        'biopython',
        'pandas',
        'seqkit',
        'bbmap'
    ]
)
