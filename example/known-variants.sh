#!/bin/bash

# Here, we simply want to create a valid subset of the 1001g VCF for testing purposes.
# We do not really care which variants make it into the final file, it just should be
# a representative subset of a size that works for exemplary and testing purposes,
# so that we can check that grenedalf and its tools generally do the right thing.
# It should include biallelic and multiallelic SNPs, as well as short indels.
# We will use that as a known variants file, as well as for training the GATK VQSR step.

# We base the file off the Arabidobis thaliana 1001 genomes VCF.
# First, we plot a histogram of the INFO DP values, to check for reasonable cutoffs.

# collect the values of the INFO DP column in a list
# zcat 1001genomes_snp-short-indel_only_ACGTN.vcf.gz | egrep -v "^#" | cut -f 8 | sed "s/DP=//g" > DP.txt

# plot a histogram for the DP values
# python3 <<END
# import matplotlib.pyplot as plt
# filename="DP.txt"
# with open(filename) as file:
#     lines = [int(line.rstrip()) for line in file]
# plt.hist(x=lines, bins="auto")
# plt.savefig('DP.png')
#
# # aux functions to do some counting of the data
# # lines.sort()
# # lines[-500000:][0]
# # sub = [ x for x in lines if x > 30000 and x < 40000 ]
# END

# List of 100 accessions from the 1001g VCF that we picked for testing here.
echo "10001,10002,10004,10005,10006,10008,10009,10010,10011,10012,10013,10014,10015,10017,10018,1002,10020,10022,10023,10027,1006,1061,1062,1063,1066,1070,108,1158,1166,1254,1257,1313,1317,139,14312,14313,14314,14315,14318,14319,1552,15560,15591,15592,15593,159,1612,1622,1651,1652,1676,1684,1739,1741,1756,1757,1793,1797,1819,1820,1829,1834,1835,1851,1852,1853,18694,18696,1872,1890,1925,1942,1943,1954,19949,19950,19951,2016,2017,2031,2053,2057,2081,2091,2106,2108,2141,2159,2166,2171,2191,2202,2212,2239,2240,2276,2278,2285,2286,2317" | sed "s/,/\n/g" > "accessions.txt"

# Now, we can subset the vcf to where the INFO field DP is somewhere between ~30k and ~40k,
# as those is the upper high quality end of the data, and yields a reasonable number
# of variants to use in our grenepipe testing.

# min and max DP INFO that we want to keep
# adjusted to give a reasonable amount of variants for testing purposes.
MIN=29000
MAX=35000

# in, out, and temp files
VCF="1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
SUB="subset.vcf"
ACC="accessions.vcf.gz"
OUT="known-variants.vcf.gz"

# select only those variants that have INFO DP between MIN and MAX
# vcftools does not seem to be able to use the INFO DP field directly for counting
# (it can use the means per sample though)... so instead we use awk :-P
zcat $VCF | head -n 100 | egrep "^#" > $SUB
zcat $VCF | egrep -v "^#" | awk -v MIN="$MIN" -v MAX="$MAX" -F'\t' 'BEGIN {OFS = FS} { DP=$8 ; sub("^DP=","",DP) ; if ( DP+0 > MIN && DP+0 < MAX ) print $0 }' >> $SUB

# select only a subset of the accessions as well
vcftools --vcf $SUB --keep "accessions.txt" --recode --recode-INFO-all --stdout | bgzip --stdout > $ACC
tabix -p vcf $ACC

# remove all annotations except for the INFO DP field, which is needed for the GATK VQSR step.
# all other tools that use this file as a known variants collection don't care about
# any annotations, so this file works for that purpose as well.
bcftools annotate -x ^INFO,^FORMAT/GT $ACC | bgzip -l 9 --stdout > $OUT
tabix -p vcf $OUT
