#!/usr/bin/env bash
# SET PATH TO REPOSITORY
REPOSITORY_PREFIX="."
#SET PATH TO CONDA ENVIRONMENT
CONDA_PREFIX=""

cd ${CONDA_PREFIX}
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
#touch ./etc/conda/activate.d/env_vars.sh
#touch ./etc/conda/deactivate.d/env_vars.sh

cp -f ${REPOSITORY_PREFIX}/envs/conda/activate/env_vars.sh ./etc/conda/activate.d/env_vars.sh
cp -f ${REPOSITORY_PREFIX}/envs/conda/deactivate/env_vars.sh ./etc/conda/deactivate.d/env_vars.sh
