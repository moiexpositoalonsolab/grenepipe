# =================================================================================================
#     FastQC
# =================================================================================================

# Make a list of all files that we want to run FastQC on.
# We store them in the global config, so that the rule can access them.
# There are several combinations of cases of which files we want to run fastqc on.
# We can have single or paired end raw read sample data; we can have the trimmed reads, which
# can also be either single or paired end, but also merged pairs. We here define a function
# that takes care of all theses cases and returns a list of all the ones we want to produce.
# The MultiQC rule then uses this list to request all those files, and the below get_fastqc_input()
# function resolves them into the input file paths.

assert "fastqc" not in config["global"]
config["global"]["fastqc"] = pd.DataFrame(columns=["sample", "unit", "id", "file"])


def add_fastqc_file(sample, unit, id, file):
    config["global"]["fastqc"] = pd.concat(
        [
            config["global"]["fastqc"],
            pd.DataFrame({"sample": [sample], "unit": [unit], "id": [id], "file": [file]}),
        ],
        ignore_index=True,
    )

    # Deprecated way of using append instead of concat. We use concat now to avoid warnings
    # with newer pandas versions, but keep the below for reference.
    # config["global"]["fastqc"] = config["global"]["fastqc"].append({
    #     "sample": sample,
    #     "unit":   unit,
    #     "id":     id,
    #     "file":   file
    # }, ignore_index=True)


if config["params"]["fastqc"]["input"] == "samples":
    # Simple case: raw fastq files from the samples.
    for smp in config["global"]["samples"].itertuples():
        add_fastqc_file(smp.sample, smp.unit, "R1", smp.fq1)
        if isinstance(smp.fq2, str):
            add_fastqc_file(smp.sample, smp.unit, "R2", smp.fq2)
elif config["params"]["fastqc"]["input"] == "trimmed":
    # Trimmed files, which can come in more varieties.
    for smp in config["global"]["samples"].itertuples():
        # We use a fake wildcard to get the function call to work.
        wc = snakemake.io.Wildcards()
        wc.sample = smp.sample
        wc.unit = smp.unit
        trimmed = get_trimmed_reads(wc)

        # Now let's see if we have merged them or not, and add to our result accordingly.
        if config["settings"]["merge-paired-end-reads"]:
            assert len(trimmed) == 1
            add_fastqc_file(smp.sample, smp.unit, "trimmed-merged", trimmed[0])
        else:
            # If not merged, it's either single end or paired end.
            assert len(trimmed) == 1 or len(trimmed) == 2
            add_fastqc_file(smp.sample, smp.unit, "trimmed-R1", trimmed[0])
            if len(trimmed) == 2:
                add_fastqc_file(smp.sample, smp.unit, "trimmed-R2", trimmed[1])
else:
    raise Exception("Unknown fastqc input setting: " + config["params"]["fastqc"]["input"])


def get_fastqc_input(wildcards):
    return (
        config["global"]["fastqc"]
        .loc[
            (config["global"]["fastqc"]["sample"] == wildcards.sample)
            & (config["global"]["fastqc"]["unit"] == wildcards.unit)
            & (config["global"]["fastqc"]["id"] == wildcards.id),
            ["file"],
        ]
        .file
    )


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="qc/fastqc/{sample}-{unit}-{id}_fastqc.html",
        zip="qc/fastqc/{sample}-{unit}-{id}_fastqc.zip",
    params:
        config["params"]["fastqc"]["extra"],
    log:
        "logs/fastqc/{sample}-{unit}-{id}.log",
    benchmark:
        "benchmarks/fastqc/{sample}-{unit}-{id}.bench.log"
    group:
        "qc"
    conda:
        "../envs/fastqc.yaml"
    script:
        # We use our own version of the wrapper here, as that wrapper is just badly implemented.
        "../scripts/fastqc.py"


# wrapper:
#     "0.27.1/bio/fastqc"


rule fastqc_collect:
    input:
        expand(
            "qc/fastqc/{u.sample}-{u.unit}-{u.id}_fastqc.html",
            u=config["global"]["fastqc"].itertuples(),
        ),
        expand(
            "qc/fastqc/{u.sample}-{u.unit}-{u.id}_fastqc.zip",
            u=config["global"]["fastqc"].itertuples(),
        ),
    output:
        touch("qc/fastqc/fastqc.done"),


localrules:
    fastqc_collect,


# =================================================================================================
#     Trimming
# =================================================================================================


# Different trimming tools produce different summary files. We here expand ourselves,
# because we need to retrieve the correct file type for each of them (single or paired end),
# which cannot easily be done with the simple snakemake expand function.
def get_trimming_reports():
    result = []
    for smp in config["global"]["samples"].itertuples():
        # The get_trimming_report() function is part of each trimming tool rule file.
        # We here hence call the respective correct function for each tool.
        result.append(get_trimming_report(smp.sample, smp.unit))

        # Now append the file for the sample to the result list
        # if config["settings"]["trimming-tool"] == "adapterremoval":
        # result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + ".settings" )
        # elif config["settings"]["trimming-tool"] == "cutadapt":
        # result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".qc-" + suffix + ".txt" )
        # elif config["settings"]["trimming-tool"] == "fastp":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + "-fastp.json" )
        # elif config["settings"]["trimming-tool"] == "skewer":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + "-trimmed.log" )
        # elif config["settings"]["trimming-tool"] == "trimmomatic":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".trimlog.log" )
        # else:
        #     raise Exception("Unknown trimming-tool: " + config["settings"]["trimming-tool"])
    return result


rule trimming_reports_collect:
    input:
        get_trimming_reports(),
    output:
        touch("trimmed/trimming-reports.done"),


localrules:
    trimming_reports_collect,
