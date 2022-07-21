#!/bin/bash

# ==================================================================================================
#      Setup
# ==================================================================================================

# Global settings
CORES=10

# Color the spectrum!
COLOR_RED="\033[31m"
COLOR_GREEN="\033[32m"
COLOR_END="\033[0m"

# Change to the parent directory of where this script is located, which is the main grenepipe
# directory, so that we can call this script from anywhere, and run everything from there.
cd `dirname ${0}`/..
BASEPATH=`pwd -P`

# Remove old output. Nope. We want to be able to continue failed tests while debugging!
# ./test/clean.sh

# Copy the reference genome and known variants from the example dir to here,
# so that we can run here without cluttering the example directory.
if [[ ! -d ./test/reference ]]; then
    mkdir -p ./test/reference/
    cp ./example/TAIR10_chr_all.fa.gz  ./test/reference/TAIR10_chr_all.fa.gz
    cp ./example/known-variants.vcf.gz ./test/reference/known-variants.vcf.gz
    cp ./example/regions.bed           ./test/reference/regions.bed
fi

# Copy the samples table, so that we can change the paths without changing the original,
# We need to use a different sed separator here that cannot occur in paths, to avoid conflict.
cat ./test/samples_template.tsv | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/samples.tsv

# For seqprep, we also want a samples table with only paired-end reads.
# So let's make a copy and remove the line that contains the single-end file.
cat ./test/samples.tsv | egrep -v "^S3" > ./test/samples-pe.tsv

# For the config, we need a bit of a special setup, as we need to replace some paths in
# different ways. We do this a bit more complicated setup here, in order to be able to work
# with the original config file, so that we do not need to keep track of copies and changes.
make_config() {
    local TARGET=$1

    cp ./config.yaml ${TARGET}
    sed -i "s?/path/to/data/samples.tsv?${BASEPATH}/test/samples.tsv?g" ${TARGET}
    sed -i "s?/path/to/data/genome.fa?${BASEPATH}/test/reference/TAIR10_chr_all.fa?g" ${TARGET}
    # cat ./test/config_template.yaml | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/config.yaml

    # Need an extra replacement step for threads. Might change in the future - this is a bit
    # volatile. But works for now.
    sed -i "s/threads: 12/threads: 6/g" ${TARGET}

    # We also copy over the VQSR training data. It's only needed for few test cases,
    # but doesn't hurt to have it in the other directories as well. As this might not be available
    # for every user that wants to run our tests, we might skip - VSQR tests fail in that case!
    if test -f "test/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"; then
        cd "$(dirname "${TARGET}")"
        ln -s ../1001genomes_snp-short-indel_only_ACGTN.vcf.gz .
        ln -s ../1001genomes_snp-short-indel_only_ACGTN.vcf.gz.tbi .
        cd - > /dev/null
    fi
}

# Helper to list _all_ decendant processes of this script, so that when killing the script,
# we can kill all of them. Killing just snakemake is not sufficient, as others keep running...
list_descendants()
{
    local children=$(ps -o pid= --ppid "$1")
    for pid in $children ; do
        list_descendants "$pid"
    done
    echo "$children"
}

# We use the snamemake version to figure out whether it supports mamba already,
# as this greatly decreases runtime when installing packages.
# The below version comparison works as long as snakemake uses a three (or four) numbering scheme,
# dot separated, without any letters or other things after the major-minor-patch versions.
CONDA_FRONTEND=""
CONDA_PREFIX="conda-envs"
function version {
    echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }'
}
if [ $(version `snakemake --version`) -ge $(version "5.18.0") ]; then
    CONDA_FRONTEND="--conda-frontend mamba"
    CONDA_PREFIX="mamba-envs"
fi

# ==================================================================================================
#      Run and Monitor Snakemake
# ==================================================================================================

# Count how many snakemake things we ran.
PASSCOUNT=0
FAILCOUNT=0

run_snakemake() {
    local CASE=$1
    local DIRECTORY=$2
    local EXTRA=$3

    # User output
    echo "[========" `date "+%F %T"` "========]"
    printf "${COLOR_GREEN}Running ${CASE}${COLOR_END}\n"

    # When we kill this script, we want snakemake (which we start in the background, see below)
    # and all its decendants to also be killed. Otherwise, it would keep running and mess up files
    # in the test dirs. Just killing the snakemake process does not work, as this keeps its children
    # alive. So, we use a hammer and kill all decendants of this script... Seems to work well so far.
    # See https://stackoverflow.com/a/8927066/4184258, https://unix.stackexchange.com/a/124148
    # and https://stackoverflow.com/a/5722874/4184258 for all sources that were needed here...
    trap 'kill -9 $(list_descendants $$) 2> /dev/null ; exit 143' TERM SIGINT SIGTERM

    # Run snakemake in the background.
    # Importantly, specify the conda prefix, so that the tools do not have to be loaded each time.
    snakemake \
        --use-conda \
        --conda-prefix ${BASEPATH}/test/${CONDA_PREFIX} \
        ${CONDA_FRONTEND} \
        --cores ${CORES} \
        --directory ${DIRECTORY} \
        ${EXTRA} \
        >> ${DIRECTORY}/test-run.log 2>&1 &
    PROC_ID=$!

    # Use the process ID to keep looping here while it is running,
    # and print some nice user output about the progress. This might not catch all progress,
    # in cases where multiple steps finish while we sleep here, but that is okay. It's for our
    # test cases only anyway, and we can live with that.
    PROGRESS=""
    while kill -0 "$PROC_ID" >/dev/null 2>&1; do
        sleep 1
        CURRENT=`egrep "[0-9]* of [0-9]* steps \([0-9.]*%\) done" ${DIRECTORY}/test-run.log | tail -n 1`
        if [[ "$CURRENT" != "$PROGRESS" ]]; then
            echo "    ${CURRENT}"
            PROGRESS=${CURRENT}
        fi
    done

    # Final user output for the test case
    SUCCESS=`grep "[0-9]* of [0-9]* steps ([0-9.]*%) done" ${DIRECTORY}/test-run.log | grep "100%"`
    NOTHING=`grep "Nothing to be done." ${DIRECTORY}/test-run.log`
    if [[ ! -z "$NOTHING" ]] ; then
        PASSCOUNT=$((PASSCOUNT+1))
        printf "${COLOR_GREEN}    Nothing to be done${COLOR_END}\n"
    elif [[ ! -z "$SUCCESS" ]]; then
        PASSCOUNT=$((PASSCOUNT+1))
        printf "${COLOR_GREEN}    Finished successfully${COLOR_END}\n"
    else
        FAILCOUNT=$((FAILCOUNT+1))
        printf "${COLOR_RED}    Error occurred${COLOR_END}\n"
        return 1
    fi
}

# ==================================================================================================
#      Test Cases
# ==================================================================================================

# Either get all scripts that we have, or the use provided ones via wildcard inclusion
DICTS=`ls ./test/cases/*.txt`
[[ "${1}" ]] && DICTS=`ls ./test/cases/*${1}*.txt`

# Now run all test cases. We always use the base config. This way, when we update the main
# grenepipe config file for new functionality, we only have to do minimal replacement here as well.
for DICT in ${DICTS} ; do

    # We collect all lines to be replaced from a "dictionary" script, one for each test,
    # and produce a new config file for this test, with all the lines replaced.
    # We also check that all the dictionary lines are actually in the original config,
    # in order to avoid that we accidentally changed our base config without noticing, which would mean
    # that we do not run the test the way we want, so better check.
    TARGET=$(basename "${DICT}" .txt)
    # echo "Running test case ${TARGET}"

    # Skip tests that start with an underscore.
    if [[ $TARGET = _* ]] ; then
        continue
    fi

    # Allow to re-run the test without doing everything again
    # (in particular, snakemake steps that already succceeded).
    # That is, do not create the directory and the config file again if those already exist.
    if [[ ! -d ./test/out-${TARGET} ]]; then

        # Make a copy of the config file with base paths correctly set up.
        mkdir ./test/out-${TARGET}
        make_config ./test/out-${TARGET}/config.yaml
        # cat ./test/config_template.yaml | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/out-${TARGET}/config.yaml

        # Do the replacement of all entries in config that are to be changed in the test case.
        while IFS="" read -r entry || [ -n "$p" ]
        do
            FROM=`echo "${entry}" | cut -f 1`
            TO=`echo "${entry}" | cut -f 2`
            # echo "from -${FROM}- to -${TO}-"

            # Need to re-echo the TO variable to resolve escape sequences.
            # Cannot do that before, as that would confuse the `cut`.
            TO=`echo -e "${TO}"`

            # Test that the FROM part exists, to make sure that we are doing the right thing
            # if we replace config lines later on.
            if ! grep -q "${FROM}" ./test/out-${TARGET}/config.yaml; then
                # For seqprep, we need special treatment: It only wants paired end sequences,
                # so we need to provide it with a different sample table. The samples table does
                # have a different path every time though, so that the grep here will fail.
                # Fetch this, and do not treat that particular one as an error.
                if [[ "${FROM}" != \samples* ]]; then
                    printf "${COLOR_RED}Test setting does not exist: ${FROM}${COLOR_END}\n"
                    exit 1
                fi
            fi

            # Replace the content in the config file. We use perl to be able to process
            # multi-line replacements as in the pileup case.
            # sed -i "s?${FROM}?${TO}?g" ./test/out-${TARGET}/config.yaml
            perl -pi -e "s?${FROM}?${TO}?g" ./test/out-${TARGET}/config.yaml
        done < ${DICT}

        # The known variants and restrict regions entries in the test cases uses a placeholder
        # for the directory, which we need to replace by the correct path here.
        sed -i "s?#BASEPATH#?${BASEPATH}?g" ./test/out-${TARGET}/config.yaml
    fi

    # For the pileup test cases, we do not want to run the whole pipeline,
    # so let's use the extra param to just run the necessary parts.
    EXTRA=""
    if [[ ${TARGET} == pileup* ]] || [[ ${TARGET} == mpileup* ]] ; then
        EXTRA="all_pileups"
    fi

    # Now run snakemake on the new config file.
    run_snakemake "${TARGET}" "./test/out-${TARGET}" "${EXTRA}"

    # Manual call of the test case. Replaced by the above wrapper call.
    # snakemake \
    #     --use-conda \
    #     --conda-prefix ${BASEPATH}/test/conda-envs \
    #     --cores ${CORES} \
    #     --directory ./test/out-${TARGET} \
    #     &> ./test/out-${TARGET}/snakemake.log
done

# Final user output
echo "[========" `date "+%F %T"` "========]"
if [[ ${PASSCOUNT} -gt 0 ]]; then
    printf "${COLOR_GREEN}PASS ${PASSCOUNT}${COLOR_END}\n"
fi
if [[ ${FAILCOUNT} -gt 0 ]]; then
    printf "${COLOR_RED}FAIL ${FAILCOUNT}${COLOR_END}\n"
fi
