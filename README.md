# README for scripts to parse barcode and insert assignments from PacBio SMRT reads

Scripts were written by Martin Kircher, the example is provided by Lea Starita.

## Circular Consensus Calling the reads

The bax2bam and ccs2 scripts are available from PacBio (https://github.com/PacificBiosciences/unanimity/). Although whomever is performing the PacBio sequencing may have done the bax2bam conversion and ccs2 circular consensus calling already. *Warning: Previous versions of the circular consensus caller do not work well enough. Be sure you are using ccs2.* 

We assume that `$PACBIO` is an environment variable with the PATH to the PacBio binaries.

You will also need samtools (http://www.htslib.org/) and some common GNU tools, which should be available in most UNIX-based environments. In our case, a version of samtools came with the PacBio scripts.

### Convert bax to bam

In our example, there are three bax.h5 files for each SMRT cell which we convert:

```
$PACBIO/bax2bam -o pacbio/NNK/bax2bam/F_bax2bam pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.1.bax.h5 pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.2.bax.h5 pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.3.bax.h5
```

### Concatenate bam files using samtools:

*Warning: BEWARE OF HEADERS*

```
$PACBIO/samtools cat -h <( ($PACBIO/samtools view -H B_bax2bam.subreads.bam; $PACBIO/samtools view -H D_bax2bam.subreads.bam; $PACBIO/samtools view -H E_bax2bam.subreads.bam; $PACBIO/samtools view -H F_bax2bam.subreads.bam ) | sort | uniq ) B_bax2bam.subreads.bam D_bax2bam.subreads.bam E_bax2bam.subreads.bam F_bax2bam.subreads.bam > bax2bam_combined.subreads.bam
```

### Run consensus calling with ccs2 

This step took days for us (despite using 32 threads). `TPMT_NNK_css2.bam` is the name of our output file.

```
$PACBIO/ccs --numThreads=32 TPMT_NNK_css2.bam bax2bam_combined.subreads.bam
```

## Aligning consensus reads to reference

We are using BWA mem for alignment. Our version was 0.7.10-r789.

Make a reference fasta  that contains all the sequence between the SMRT bell adaptors, *plus a few bases* of adaptor. *Delete the barcode.* The barcode sequence will then become an insertion in the BWA alignment

Example:

```
vim TPMT.fa
```

```
>TPMT_nobarcode
TGTACAAGATGGATGGTACAAGAACTTCACTTGACATTGAAGAGTACTCGGATACTGAGGTACAGAAAAACCAAGTACTAACTCTGGAAGAATGGCAAGACAAGTGGGTGAACGGCAAGACTGCTTTTCATCAGGAACAAGGACATCAGCTATTAAAGAAGCATTTAGATACTTTCCTTAAAGGCAAGAGTGGACTGAGGGTATTTTTTCCTCTTTGCGGAAAAGCGGTTGAGATGAAATGGTTTGCAGACCGGGGACACAGTGTAGTTGGTGTGGAAATCAGTGAACTTGGGATACAAGAATTTTTTACAGAGCAGAATCTTTCTTACTCAGAAGAACCAATCACCGAAATTCCTGGAACCAAAGTATTTAAGAGTTCTTCGGGGAACATTTCATTGTACTGTTGCAGTATTTTTGATCTTCCCAGGACAAATATTGGCAAATTTGACATGATTTGGGATAGAGGAGCATTAGTTGCCATTAATCCAGGTGATCGCAAATGCTATGCAGATACAATGTTTTCCCTCCTGGGAAAGAAGTTTCAGTATCTCCTGTGTGTTCTTTCTTATGATCCAACTAAACATCCAGGTCCACCATTTTATGTTCCACATGCTGAAATTGAAAGGTTGTTTGGTAAAATATGCAATATACGTTGTCTTGAGAAGGTTGATGCTTTTGAAGAACGACATAAAAGTTGGGGAATTGACTGTCTTTTTGAAAAGTTATATCTACTTACAGAAAAGTAACTCGAGCATATGACATGTCCTAGGCTTAAGCTAGCTCTAGACTGATCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGC
```

### Create index file for reference fasta

```
$PACBIO/samtools faidx ../reference/TMPT.fa
```

### Create reference genome using reference fasta with BWA:

```
bwa index -a is ../reference/TMPT.fa
```

### Run BWA alignment:

Example shell script:

```
#!/bin/bash
bwa mem -C -M -L 80 ../reference/TMPT.fa <($PACBIO/samtools view -F 1 ../TPMT_NNK_css2.bam | awk 'BEGIN{ FS="\t"; OFS="\n" }{ split($0,a,"\t"); helper = ""; for (i=12; i <= length(a); i++) { helper = helper""a[i]"\t"}; sub("\t$","",helper); print "@"$1" "helper,$10,"+",$11 }' ) | $PACBIO/samtools view -uS - | $PACBIO/samtools sort -o TPMT_NNK_css2_aligned.bam - 
$PACBIO/samtools flagstat TPMT_NNK_css2_aligned.bam > TPMT_NNK_css2_aligned.bam_stats
```

### Check alignment CIGARS to make sure your reference is correct

Check CIGARS to make sure reference is correct and there was enough extra sequence on the ends of the reference to avoid soft-clipping. There should be a majority CIGAR with matches to the whole sequence and an insertion the size of the barcode. Very important not to have a majorty of insertions or soft-clips on either end!

```
$PACBIO/samtools view TPMT_NNK_css2_aligned.bam | awk '{ if ($4 == 1) print $6 }' | sort | uniq -c | sort -nr | head
44091 781M15I221M      ###Read as: 781 matches, 15 base insertion, 221 matches
13928 780M15I222M
3208 781M15I162M1D58M
2518 779M15I223M
990 780M15I163M1D58M
875 778M15I224M
792 781M15I189M1D31M
655 302M1D478M15I221M
648 1002M
619 781M14I221M
```

### Remove PacBio specific columns from BAM file

Cutoff a bunch of PacBio specific columns from the bam file so it is compatible with pysam (everything after field 13):

```
$PACBIO/samtools view -h TPMT_NNK_css2_aligned.bam.bam | cut -f -14 | $PACBIO/samtools view -Sb - > TPMT_NNK_css2_aligned.fix.bam
```

## Parsing alignments to identify barcode and associated insert sequence

We used the following software versions:
python/2.7.3
pysam/0.7.5

`extractBarcodeInsertPairs.py` parses the CIGAR and MD strings to identify the insert (variable region) and the barcode. See the parameters below and make sure to adjust the defaults (e.g. `length`, `position`, `start`, and `end`) for your specific project. The -v option will let you know which filters sequences are passing.

```
python extractBarcodeInsertPairs_moreQC.py --help

Options:
  -h, --help            show this help message and exit
  -l BLENGTH, --length=BLENGTH
                        Length of the barcodes (default 16)
  -p BPOS, --position=BPOS
                        Position of the barcode (default 52)
  --start=START         Start reference position for extracted sequence
                        (default 114)
  --end=END             End reference position for extracted sequence (default
                        1065)
  -b MINBASEQ, --minBaseQ=MINBASEQ
                        Minimum base quality score (default 0)
  -q MINMAPQ, --minMapQ=MINMAPQ
                        Filter alignments by mapQ (default 0)
  -s, --allow-soft-clipped
                        Consider soft clipped reads (default Off)
  -v, --verbose         Turn debug output on
```

In our example:

```
python extractBarcodeInsertPairs_moreQC.py TPMT_NNK_css2_aligned.fix.bam -l 15 -p 781 --start=8 --end=746 | gzip -c > TPMT_NNK_barcodeInsertQuals.tsv.gz
```

## Unifying barcode reads

unifyAssignment.py assigns the most frequently associated sequence or the sequence with the best average quality score to each barcode.

Command line used in our example:

```
python unifyAssignment.py TPMT_NNK_barcodeInsertQuals.tsv.gz | gzip -c > TMPT_NNK_barcodeInsertAssignment.tsv.gz
```
