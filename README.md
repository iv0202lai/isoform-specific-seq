# Isoform-specific-seq
by Ivan Pochou Lai at School of Medicine, National Taiwan University

## Introduction
This is a tool used to get isoform-specific sequence of mice from an Ensembl transcript number, that is, to find out the sequence unique to other isoforms of the same gene.
The tool is also capable of marking SNP sites for B6/CAST on the output sequence in the mean while.

Although the project was meant to mark the SNP sites of B6/CAST in mice, other species and SNP reference files should work with a minor modification.

## Prerequisites
* Linux or MacOS
* Python 2
* PyCogent 1.9

## Input
1. Simply enter an Ensembl transcrpit number.
2. Use a txt file that contains multiple Ensembl transcrpit numbers, where each line contains a single transcript number.

## Output
1. Shown on the screen.
2. Saved as txt files in a folder called "isoform_specific_seq" under local folder.

## SNP site marking for B6/CAST comparison
You need 2 files named "CAST_snps.vcf.gz" and "CAST_snps.vcf.gz.tbi" as reference files.
You can use your own files or download the ones from the [website of Mouse Genome Project of Sanger's institute](https://www.sanger.ac.uk/science/data/mouse-genomes-project) .
Put the referecne file in the same folder and SNP sites will be automatically marked on the result sequence in the form of \[B6/CAST\].
If there're no designated reference files, the program will show the output without SNP site labeling.

Sample format of vcf files (extracted from "CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf" on Sanger's instute's FTP)
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CAST_EiJ
1	3000023	.	C	A	39.2026	MinDP;MinAB	DP4=0,0,4,0;DP=4;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:15:4:0:79,15,0:67,12,0:2:24:4:0,0,4,0:0:-0.556411:.:0
1	3000185	rs585444580	G	T	228	PASS	DP4=0,0,10,2;DP=12;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:43:12:0:276,43,0:255,36,0:2:54:12:0,0,10,2:0:-0.680642:.:1
```

## Example
Input an Ensembl mouse transcrpit number
```
ENSMUST00000079360
```
The output .txt file will show
```
Transcript ID: ENSMUST00000079360
Gene: Ablim1(ENSMUSG00000025085)
Chromosome: 19       Strand: Reverse

Region 1   Position: 57032733-57032889   Length: 157   SNP count: 4
[C/T]CA[G/C]TGCTGAGCACCTGGTATTG[A/G]GGTACAGAGGAGGGCTCCGCCTGCAGGCACTGTGGCCGCTGAGCTGCGTGTTCCCAGCTGCCTCTGTTGCGGACAGCCCGCT[C/T]CACTCTCCACCCTGTACTCCATAAATAAAGCTACGTGCTTCTGCCCTGAG

Region 2   Position: 57215966-57216032   Length: 67   no SNP
AGGGGAACTGCTACTGGGGAGCCCTGCTTCCCGGGCTTCTGGGTCAGATCTTGATGTCATGTCTGCT
```
You can also input a txt file contain a set of transcript numbers
```
ENSMUST00000154958
ENSMUST00000144949
ENSMUST00000079360
```
There will be 3 output .txt files named as their transcript numbers in ./isoform_specific_seq
