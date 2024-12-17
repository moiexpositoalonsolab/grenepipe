Grenepipe External Dependencies
=============================

This directory contains packages (external dependencies) of grenepipe
that are not available via the conda package manager (i.e., are not on conda-forge or bioconda).
We currently include

 - [HAFpipe](https://github.com/lczech/HAF-pipe), to be downloaded into `hafpipe`. Note that we use our own form, which contains several improvements over the original [HAFpipe-line](https://github.com/petrov-lab/HAFpipe-line).
 - [harp](https://bitbucket.org/dkessner/harp/), to be downloaded into `harp`. Here, we only need the binary, and do not need the repository to be cloned.

Calling

    ./setup-hafpipe.sh

in this directory will download the dependencies and make them available for grenepipe.
However, this is also done automatically when running grenepipe with a config file that
requests HAFpipe to be run; we want life to be easy for you!
