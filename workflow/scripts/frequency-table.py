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

# Get the columns that we want to print. We do some extra work to ensure the exact order that we want.
# That's needed so that our printing below prints out the correct order without extra work there.
fields_tmp = snakemake.params.get("fields", "COV,FREQ,REF_CNT,ALT_CNT").split(",")
fields = []
for f in ["COV", "FREQ", "REF_CNT", "ALT_CNT"]:
    f = f.strip()
    if f in fields_tmp:
        fields.append(f)
        fields_tmp.remove(f)
if len(fields_tmp) > 0:
    print("Invalid fields are ignored: " + str(fields_tmp))

# Make header with all combinations of these fields for all samples.
# We tried itertools.permutations, but that expects the longer list to go first and all that.
# Nope, doing it by hand then... Maybe I should learn proper python at some point...
sample_fields = []
for field in fields:
    for sample in reader.header.samples.names:
        sample_fields.append(field + "." + sample)

# Print header
header = ["CHROM", "POS", "REF", "ALT"] + sample_fields
tableout.write("\t".join(header) + "\n")

# Process input vcf line by line
for record in reader:
    if not record.is_snv() or len(record.ALT) > 1:
        continue
    # print(str(record))

    # Prepare fixed part
    line = [record.CHROM, record.POS, record.REF]
    line += [alt.value for alt in record.ALT]
    tableout.write("\t".join(map(str, line)))

    # Prepare lists for results per sample
    covs = []
    freqs = []
    ref_cnts = []
    alt_cnts = []

    # Iterate samples and fill lists
    # line += [call.data.get('AD') or '0,0' for call in record.calls]
    for call in record.calls:
        ad = call.data.get("AD")
        if ad and len(ad) == 2:
            cov = ad[0] + ad[1]
            covs.append(str(cov))
            ref_cnts.append(str(ad[0]))
            alt_cnts.append(str(ad[1]))
            if cov > 0:
                frac = ad[0] / cov
                freqs.append(str(frac))
            else:
                freqs.append("NA")
        else:
            # Fill with NA if there is no AD in the VCF
            covs.append("NA")
            freqs.append("NA")
            ref_cnts.append("NA")
            alt_cnts.append("NA")

    # Write the lists per sample, in the order of the fields.
    if "COV" in fields:
        tableout.write("\t" + "\t".join(covs))
    if "FREQ" in fields:
        tableout.write("\t" + "\t".join(freqs))
    if "REF_CNT" in fields:
        tableout.write("\t" + "\t".join(ref_cnts))
    if "ALT_CNT" in fields:
        tableout.write("\t" + "\t".join(alt_cnts))
    tableout.write("\n")
tableout.close()
