# MADT
**Meta Amplicon Database Toolkit (MADT)** attempts to provide a standardized preprocessing process for the current massive open-source amplification data. In the case of missing or inaccurate meta-information, it screens for errors in common sequence files and retains the correct sequence files, ultimately establishing a standardized in-slico DNA library. The standardized in-slico DNA library can largely benefit downstream meta-analysis. Currently, MADT is still in the internal testing phase. If you want to experience basic functions or participate in MADT testing, please contact yfzeng0827@hotmail.com.
## Install by pip
Now MADT only can be installed by pip:
```
pip install
```
MADT only can be used in Linux command line. Please install *bbmaps* and *seqkit* for the dependency.

## Brief tutorial
(We will give more details in the future version)
(i) We need to prepare a file with SRA accession ID which we want to download line by line.
```
madt_DataIndex.py -i SRA.list.txt -o raw_data
```
(ii) We sort out the potential single-end files or pair-end files. In this step, low-quality files will be discarded. Different files will be stored in different output directory
```
madt_DataIndex.py -i raw_data
```
(iii) We do a quality-control step on single-end files or pair-end files. Pair-end files will be merged. 
```
madt_LibCheck.py -s single_data -p pair_data
```
(iV) We use blast to determine the amplicon region of files. And sort out files by amplicon region. Now only support to detect bacteria. Please use the reference file in the *blast* directory of our packages.
```
madt_LibMap.py -i pair_data -b blast -r Ref/Ecoil/Ecoil
madt_LibMap.py -i single_data -b blast -r Ref/Ecoil/Ecoil
```