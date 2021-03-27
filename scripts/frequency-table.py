# =================================================================================================
#     Dependencies
# =================================================================================================

import vcfpy

# =================================================================================================
#     Process Data
# =================================================================================================

# Prepare output file
tablefile = snakemake.output[0]
tableout = open(tablefile, "w")

# Prepare input file
vcffile = snakemake.input[0]
reader = vcfpy.Reader.from_path(vcffile)

# Print header
header = ['CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names
tableout.write('\t'.join(header) + '\n')

# Process input vcf line by line
for record in reader:
    if not record.is_snv() or len(record.ALT) > 1:
        continue
    # print(str(record))

    # Prepare fixed part
    line = [record.CHROM, record.POS, record.REF]
    line += [alt.value for alt in record.ALT]
    tableout.write('\t'.join(map(str, line)))

    # Iterate samples
    # line += [call.data.get('AD') or '0,0' for call in record.calls]
    for call in record.calls:
        ad = call.data.get('AD')
        if ad and len(ad) == 2:
            sum = ad[0] + ad[1]
            if sum > 0:
                frac = ad[0] / sum
                tableout.write('\t' + str(frac))
            else:
                tableout.write('\tNA')
        else:
            tableout.write('\tNA')

    tableout.write('\n')
tableout.close()
