This is a tool used to get isoform-specific sequence of mice from an Ensembl transcript number, that is, to find out the sequence unique to other isoforms of the same gene.
The tool is also capable of marking SNP sites for B6/CAST on the output sequence in the mean while.
Codes are written in python 2 and use PyCogent 1.9 to connect to Ensembl database

# Usage
## Input type
1. Simply enter an Ensembl transcrpit number
2. Use a txt file that contains multiple Ensembl transcrpit numbers, where each line contains a single transcript number

## Output type
1. Shown on the screen.
2. Out put as txt files. The files will be saved in a folder called "isoform_specific_seq" under local folder.

## SNP site marking for B6/CAST comparison
You need 2 files named "CAST_snps.vcf.gz" and "CAST_snps.vcf.gz.tbi" as reference file, which can be downloaded from the [website of Mouse Genome Project of Sanger's institute](https://www.sanger.ac.uk/science/data/mouse-genomes-project) 
Put the referecne file in the same folder and SNP sites will be automatically marked on the result sequence in the form of \[B6/CAST\].
If there're no designated reference file, the program will show the output without SNP site labeling.

Although the project was meant to mark the SNP sites of B6/CAST in mice, other species and SNP reference files should work with a minor modification.

# Example
Given the input
```
ENSMUST00000079360
```
The output will be
```
Transcript ID: ENSMUST00000079360
Gene: Ablim1(ENSMUSG00000025085)
Chromosome: 19       Strand: Reverse

Region 1   Position: 57032733-57032889   Length: 157   SNP count: 4
[C/T]CA[G/C]TGCTGAGCACCTGGTATTG[A/G]GGTACAGAGGAGGGCTCCGCCTGCAGGCACTGTGGCCGCTGAGCTGCGTGTTCCCAGCTGCCTCTGTTGCGGACAGCCCGCT[C/T]CACTCTCCACCCTGTACTCCATAAATAAAGCTACGTGCTTCTGCCCTGAG

Region 2   Position: 57215966-57216032   Length: 67   no SNP
AGGGGAACTGCTACTGGGGAGCCCTGCTTCCCGGGCTTCTGGGTCAGATCTTGATGTCATGTCTGCT
```
