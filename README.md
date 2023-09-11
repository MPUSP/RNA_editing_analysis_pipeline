# RNA editing analysis pipeline for SE and PE NGS data
RNA editing pipeline / code for publication Wulff et al.
Please cite the following publication when using this pipeline:<br>

## Description:
This pipeline is designed to analyse RNA editing in paired-end sequencing data and it is designed to work with DNA and RNA sequencing data from the same sample (paired DNA and RNA data). Paired-end sequencing data for input is recommended, but is optional and the pipeline can run with forward data only and with mixed datasets.

## Important:
This code has been developed to conduct analyses as outlined in the publication above and is currently in development. The code might contain bugs and is also not written for speed or for readability. Future releases might improve on those aspects.

## Disclaimer:

The software "RNA_editing_analysis_pipeline" is intended for research purposes only and is provided "as is" without warranty of any kind. The developer(s) and provider(s) of this software make no representations or warranties, express or implied, regarding the use or performance of this software.

You assume full responsibility and risk for the use of "RNA_editing_analysis_pipeline". The developer(s) and provider(s) shall not be liable for any direct, indirect, incidental, special, exemplary, or consequential damages arising out of the use of this software.

"RNA_editing_analysis_pipeline" is not intended to be used for critical, medical, or life-saving purposes. It should not be relied upon as the sole basis for decision-making.

By using "RNA_editing_analysis_pipeline", you acknowledge and agree to this disclaimer. If you do not agree to these terms, you should refrain from using the software.

Last updated: 2023-09-07

## Copyright:
Code has been created by Knut Finstermeier, Max Planck Unit for the Science of Pathogens, Berlin, Germany, member of the Max Planck Society, https://www.mpg.de/en. The code is licensed under the MIT license (see LICENSE file).

## Release:
Version 6.0, 1.8.2023

## Pre-requisites (versions are suggestions and others may work):
- Python 3.9 (https://www.python.org/downloads/release/python-390/)
- incompatible with Python 3.11 due to libffi-dev being required by pysam and cPython
- UMI-tools 1.1.0 (https://github.com/CGATOxford/UMI-tools)
- STAR aligner 2.6.1d (https://github.com/alexdobin/STAR)
- BWA aligner 0.7.17 (https://github.com/lh3/bwa)
- SAMtools 1.9 (https://github.com/samtools/samtools)
- Cutadapt 2.10 (https://github.com/marcelm/cutadapt)

Additionally, the following reference fasta sequence is required:
- PhiX (https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1)

## Installation:
1. Clone repository<br>```git clone https://github.com/MPUSP/RNA_editing_analysis_pipelin.git```
2. Run interactive ```Installer.sh``` with user permissions from within the repository. The installer will ask for the paths to the requirements above (use tab for path completion).
3. Run the pipeline script just with the parameter -h to see the help message and the required input parameters or read the description below. This will further confirm that all Python libraries are installed correctly.

## To run:
1. In CLI, create your output directory, i.e. ```mkdir test_run``` and enter into that directory ```cd test_run```
2. Provide the input file list (4 columns, tab-separated) in the following format:<br>
```DNA or RNA | matched sample name | fwd or rev | full path to fastq file```<br>
Example:<br>
```DNA    1A	fwd	/home/user/1A_fwd.fastq.gz```<br>
```DNA    1A	rev	/home/user/1A_rev.fastq.gz```<br>
```RNA    1A	fwd	/home/user/1A_fwd.RNA.fastq.gz```<br>
```RNA    1A	rev	/home/user/1A_rev.RNA.fastq.gz```<br>

3. Run the script from within your location by using i.e. absolute paths ```/path/to/script/rna_editing_full_analysis.06.py```
4. As a minimum, the script requires the path to the input file list ```-f```, the path to the PhiX fasta file ```-x```, the path to the reference fasta file ```-r``` and the reference annotation file (gff format) ```-a```.

## Generic run example:
```./rna_editing_full_analysis.06.py -f input_file_list.txt -x references/PhiX.fasta -r S_pyogenes.fasta -a S_pyogenes_anno.gff -n 10 -o test_run```<br>
This will run the pipeline with ```10``` CPUs, processing the sample read data from ```input_file_list.txt```, use ```S_pyogenes.fasta``` as fasta reference sequence, ```S_pyogenes_anno.gff``` as reference annotation, map reads against the PhiX sequence, and will create the output with the prefix ```test_run```.

## All parameters:
#### Regular parameters:
| Parameter | Description [default] |
| --- | --- |
| -h | help |
| --version | returns current version |
| -f | tab file with list of sequencing read input files [""] |
| -r | reference fasta [""] |
| -a | reference annotation file, gff3 format [""]<br> |
| -x | PhiX fasta sequence [""] |
| -n | number of CPUs for parallel processing [1] |
| -o | output filename prefix, can contain paths [""] |
| --cov | minimun coverage for called sites [5] |
| --direction | minimum reads per orientation for support of sites [1] |

#### Development parameters:
| Parameter | Description [default] |
| --- | --- |
| -d | developer mode (sets some parameters used during development) [OFF] |
| --nomap | skip mapping step if data structure exists already [False] |
| -l | log file, only partially implemented ["pipeline.log"] |
| -s | omits analysis steps (binary format, i.e. 000111) [all on] |
| --att | creates read attrition plot [on] |
| --no_umi | omits UMIs in reads (i.e. when absent) and uses regular read deduplication |
| --term | terminates analysis prematurely after this provided step [0=none] |
| --askip | PARTIALLY IMPLEMENTED: skips existing and complete steps in analysis folder |

## Acknowledgements:
The inspiration for the code came from the following publication:<br>
Bar-Yaacov, Dan et al. “RNA editing in bacteria recodes multiple proteins and regulates an evolutionarily conserved toxin-antitoxin system.” Genome research vol. 27,10 (2017): 1696-1703. doi:10.1101/gr.222760.117
Labwork has been conducted by Karin Hahnke, and Thomas F. Wulff. Sequencing was performed by the Max Planck Sequencing Core Unit at the MPI for Molecular Genetics in Berlin, Germany.

## Contact:
Please use software-dev@mpusp.mpg.de for any questions or comments.
