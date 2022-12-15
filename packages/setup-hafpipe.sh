#!/usr/bin/env bash

# Change to the directory where this script is located,
# so that we are able to call it from anywhere.
cd `dirname ${0}`

# HAFpipe is included from GitHub, and harp is directly downloaded as a binary,
# as building the binary from its sources is near impossible due to its old code base
# and un-updated C++ code that just fails to compile.
#
# We also played around with git submodules:
#
#     git submodule add https://github.com/petrov-lab/HAFpipe-line.git hafpipe
#     git submodule add https://bitbucket.org/dkessner/harp harp
#
# We removed harp as a submodule again, because we need its binary anyway,
# and then did the same for HAFpipe, to keep it simple and not have to use submodules.

# -------------------------------------------------------------------------
#     Download harp
# -------------------------------------------------------------------------

# If the harp binary is not there, download it.
HARP_BIN="harp/bin/harp"
if [ ! -f "${HARP_BIN}" ] ; then

    echo "=================================================="
    echo "    Downloading harp"
    echo "=================================================="
    echo

    # Download the respective binary of harp for the OS we are on.
    # Switching between the OS'es: https://stackoverflow.com/a/17072017/4184258
    if [ "$(uname)" == "Darwin" ]; then
        # Mac OS X
        HARP_ZIP="harp_osx_140925_100959.zip"

    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        # GNU/Linux
        HARP_ZIP="harp_linux_140925_103521.zip"

    else
        # Neiter Mac nor Linux
        echo "harp binaries are only available for Mac and Linux."
        echo "Please visit https://bitbucket.org/dkessner/harp to build your own binary."
        echo "Then, store it in `pwd`/harp/bin for grenepipe to be able to find it."
        return 1
    fi

    # We use curl instead of wget here, as this seems to be available on both Mac and Linux.
    # We need to tell curl to follow redictions... see https://askubuntu.com/a/1036492/561263
    if [ ! -f "${HARP_ZIP}" ] ; then
        curl -LO "https://bitbucket.org/dkessner/harp/downloads/${HARP_ZIP}" -o "${HARP_ZIP}"
    fi
    if [ ! -d "${HARP_ZIP%.zip}" ] ; then
        unzip "${HARP_ZIP}"
    fi

    # Unfortunately, the zip files already contain a directory called after the filename itself...
    # so either "harp_linux_140925_103521" or "harp_osx_140925_100959"; let's rename this,
    # so that we always end up having harp in the "harp" dir.
    # Also, when this script is called from the snakemake rule, snakemake is so nice as to create
    # all directories for us already, so that the below intended renaming via mv instead moves
    # the directory into the snakemake-create one... So we need to delete that one first.
    if [ -d "harp" ] ; then
        rm -r "harp"
    fi
    mv "${HARP_ZIP%.zip}" "harp"

    # If curl breaks down in the future because of Atlassion (bitbucked) changing their API,
    # try the following, which is adapted from
    # https://community.atlassian.com/t5/Bitbucket-questions/curl-command-downloading-false-zip-file-from-bitbucket-after-1/qaq-p/1973932
    # curl -X GET -O -L \
    #   --url https://api.bitbucket.org/2.0/repositories/dkessner/harp/downloads/${HARP_ZIP} \
    #   --header 'Authorization: Basic AppPasswordHere'

    if [ ! -f "${HARP_BIN}" ] ; then
        echo "Could not download harp."
        echo "Please visit https://bitbucket.org/dkessner/harp to download manually, "
        echo "using the OSX or Linux binary from 2014-09-25."
        echo "Then, store the binary in `pwd`/harp/bin for grenepipe to be able to find it."
        return 1
    fi
    echo
fi

# -------------------------------------------------------------------------
#     Download HAFpipe
# -------------------------------------------------------------------------

# If the HAFpipe binary is not there, download it.
HAFPIPE_BIN="hafpipe/HAFpipe_wrapper.sh"
if [ ! -f "${HAFPIPE_BIN}" ] ; then

    echo "=================================================="
    echo "    Downloading HAFpipe"
    echo "=================================================="
    echo

    # We use the (as of now) latest commit of HAFpipe, but still want to use a fixed commit,
    # to ensure future compatibility. We can manually update if needed.

    # Original
    # HAFPIPE_USER="petrov-lab"
    # HAFPIPE_REPO="HAFpipe-line"
    # HAFPIPE_HASH="0b773910ee625aacc5c7a43f9cf1a7f0eb6c5da2"

    # Our improved version
    HAFPIPE_USER="lczech"
    HAFPIPE_REPO="HAF-pipe"
    HAFPIPE_HASH="ef09729327bd305bdc36cb9e34462a9c85da4f75"

    # Somehow, curl from GitHub works differently - the output is already put into a file,
    # and we get a warning "Warning: Got more output options than URLs" when specifying `-o`,
    # so for now, we don't. If this does not work for some users, we'll have to revisit this.
    # See also https://gist.github.com/jwebcat/5122366 for the `J` speficier.
    HAFPIPE_ZIP="${HAFPIPE_REPO}-${HAFPIPE_HASH}.zip"
    if [ ! -f "${HAFPIPE_ZIP}" ] ; then
        curl -LJO "https://github.com/${HAFPIPE_USER}/${HAFPIPE_REPO}/archive/${HAFPIPE_HASH}.zip"
    fi
    # curl -LO "https://github.com/${HAFPIPE_USER}/${HAFPIPE_REPO}/archive/${HAFPIPE_HASH}.zip" -o "${HAFPIPE_HASH}.zip" -

    # We use curl with -J, which uses the server-side name for the downloaded file,
    # which differs from the filename as given in the URL above, but is the actual name of the file.
    # This ensures that the unzip command does not complain about mismatching file names
    # between the zip file name and the directory contained in it.
    if [ ! -d "${HAFPIPE_ZIP%.zip}" ] ; then
        unzip "${HAFPIPE_ZIP}"
    fi

    # We want to rename this to a path that we can use in grenepipe.
    if [ -d "hafpipe" ] ; then
        rm -r "hafpipe"
    fi
    mv "${HAFPIPE_REPO}-${HAFPIPE_HASH}" "hafpipe"

    if [ ! -f "${HAFPIPE_BIN}" ] ; then
        echo "Could not download HAF-pipe."
        echo "Please visit https://github.com/${HAFPIPE_USER}/${HAFPIPE_REPO} to download manually, "
        echo "using the version at commit hash ${HAFPIPE_HASH}."
        echo "Then, store the files in `pwd`/hafpipe for grenepipe to be able to find it."
        return 1
    fi
    echo
fi
