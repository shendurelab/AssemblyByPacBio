##########################################
#### README for using Martin Kircher's ###
#### scripts to parse barcode and      ###
#### insert from PacBio SMRT reads     ###
##########################################


##########################################
### Circular Consensus Calling the reads##
##########################################

The bax2bam and ccs2 scripts are available from PacBio  https://github.com/PacificBiosciences/unanimity/
Although whomever is performing the PacBio sequencing may be able to do the bax2bam conversion and ccs2 circular consensus calling. 
## warning the previous versions of the circular consensus caller do not work well enough to use for this purpose be sure you are using ccs2. 

convert bax to bam
There are three bax.h5 files for each SMRT cell
example from our work:

$PACBIO/bax2bam -o /net/shendure/vol10/projects/GPS/data/pacbio/NNK/bax2bam/F_bax2bam /net/shendure/vol10/projects/GPS/data/pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.1.bax.h5 /net/shendure/vol10/projects/GPS/data/pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.2.bax.h5 /net/shendure/vol10/projects/GPS/data/pacbio/NNK/F01_1/Analysis_Results/m160713_210035_42134_c101034302550000001823250511171615_s1_p0.3.bax.h5 &

concatenate bam files using samtools, BEWARE OF HEADERS:

example:
$PACBIO/samtools cat -h <( ($PACBIO/samtools view -H B_bax2bam.subreads.bam; $PACBIO/samtools view -H D_bax2bam.subreads.bam; $PACBIO/samtools view -H E_bax2bam.subreads.bam; $PACBIO/samtools view -H F_bax2bam.subreads.bam ) | sort | uniq ) B_bax2bam.subreads.bam D_bax2bam.subreads.bam E_bax2bam.subreads.bam F_bax2bam.subreads.bam > bax2bam_combined.subreads.bam


run ccs2  (it took days)
(TPMT_NNK_css2.bam = outfile)

$PACBIO/ccs --numThreads=32 TPMT_NNK_css2.bam bax2bam_combined.subreads.bam

############################################
### aligning consensus reads to reference ##
############################################
software versions:
bwa version: 0.7.10-r789
samtools version: 0.1.18-evan.3

Make a reference fasta  that contains all the sequence between the SMRT bell adaptors, plus a few bases of adaptor. Delete the barcode. The barcode sequence will then become an insertion in the BWA alignment

example:

vim TPMT.fa

>TPMT_nobarcode
TGTACAAGATGGATGGTACAAGAACTTCACTTGACATTGAAGAGTACTCGGATACTGAGGTACAGAAAAACCAAGTACTAACTCTGGAAGAATGGCAAGACAAGTGGGTGAACGGCAAGACTGCTTTTCATCAGGAACAAGGACATCAGCTATTAAAGAAGCATTTAGATACTTTCCTTAAAGGCAAGAGTGGACTGAGGGTATTTTTTCCTCTTTGCGGAAAAGCGGTTGAGATGAAATGGTTTGCAGACCGGGGACACAGTGTAGTTGGTGTGGAAATCAGTGAACTTGGGATACAAGAATTTTTTACAGAGCAGAATCTTTCTTACTCAGAAGAACCAATCACCGAAATTCCTGGAACCAAAGTATTTAAGAGTTCTTCGGGGAACATTTCATTGTACTGTTGCAGTATTTTTGATCTTCCCAGGACAAATATTGGCAAATTTGACATGATTTGGGATAGAGGAGCATTAGTTGCCATTAATCCAGGTGATCGCAAATGCTATGCAGATACAATGTTTTCCCTCCTGGGAAAGAAGTTTCAGTATCTCCTGTGTGTTCTTTCTTATGATCCAACTAAACATCCAGGTCCACCATTTTATGTTCCACATGCTGAAATTGAAAGGTTGTTTGGTAAAATATGCAATATACGTTGTCTTGAGAAGGTTGATGCTTTTGAAGAACGACATAAAAGTTGGGGAATTGACTGTCTTTTTGAAAAGTTATATCTACTTACAGAAAAGTAACTCGAGCATATGACATGTCCTAGGCTTAAGCTAGCTCTAGACTGATCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGC

Create index file for reference fasta:
example:
$PACBIO/samtools faidx ../reference/TMPT.fa

Create reference genome using reference fasta with BWA:
example:
bwa index -a is ../reference/TMPT.fa

Run BWA alignment:
example shell script:

#!/bin/bash                                                                                                                                                                                                                                  
bwa mem -C -M -L 80 ../reference/TMPT.fa <(samtools view -F p ../TPMT_NNK_css2.bam | awk 'BEGIN{ FS="\t"; OFS="\n" }{ split($0,a,"\t"); helper = ""; for (i=12; i <= length(a); i++) { helper = helper""a[i]"\t"}; sub("\t$","",helper); print "@"$1" "helper,$10,"+",$11 }' ) | samtools view -uS - | samtools sort - TPMT_NNK_css2_aligned.bam 
samtools flagstat TPMT_NNK_css2_aligned.bam > TPMT_NNK_css2_aligned.bam_stats

Check CIGARS to make sure reference is correct and there was enough extra sequence on the ends of the reference to avoid soft-clipping. There should be a majority CIGAR with matches to the whole sequence and an insertion the size of the barcode. Very important not to have insertions or soft-clips on either end!
example:
samtools view TPMT_NNK_css2_aligned.bam | awk '{ if ($4 == 1) print $6 }' | sort | uniq -c | sort -nr | head
44091 781M15I221M      ###781 matches, 15 base insertion, 221 matches
13928 780M15I222M
3208 781M15I162M1D58M
2518 779M15I223M
990 780M15I163M1D58M
875 778M15I224M
792 781M15I189M1D31M
655 302M1D478M15I221M
648 1002M
619 781M14I221M

Cutoff a bunch of PacBio specific columns from the bam file so it is compatible with pysam (everything after field 13):
example:
$PACBIO/samtools view -h TPMT_NNK_css2_aligned.bam.bam | cut -f -14 | $PACBIO/samtools view -Sb - > TPMT_NNK_css2_aligned.fix.bam


#############################################
### parsing alignments to identify barcode ##
###            and insert 				   ##
#############################################
software versions:
python/2.7.3
pysam/0.7.5
gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1

Martin Kircher wrote extractBarcodeInsertPairs to parse the CIGAR and MD strings to identify insert (variable region) and the barcode. Parameters below. The -v option will let you know which filters sequences are passing.

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


example:
python extractBarcodeInsertPairs_moreQC.py.py TPMT_NNK_css2_aligned.fix.bam -l 15 -p 781 --start=8 --end=746 | gzip -c > TPMT_NNK_barcodeInsertQuals.tsv.gz


#############################################
###   unifying barcode reads			   ##
#############################################

Martin Kircher wrote unifyAssignment.py to assign the most associated sequence or sequence with the best average quality score to each barcode

example:
python unifyAssignment.py TPMT_NNK_barcodeInsertQuals.tsv.gz | gzip -c > TMPT_NNK_barcodeInsertAssignment.tsv.gz



