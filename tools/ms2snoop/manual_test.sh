#!/bin/bash

VENV=.venv
CENV=.conda
MAMBA_SOLVER=
## comment to deactivate mamba solver
MAMBA_SOLVER="--experimental-solver libmamba"

SIRIUS_VERSION=4.9.15

ENV_NAME=ms2snoop_2_0_0
XML_PATH=./MS2snoop.xml

if [ ! -e "./${VENV}" ];then
  echo "virtualenv not created yet, creating..."
  python3 -m virtualenv "${VENV}"
  echo "venv created"
else
  echo "virtualenv already exist: ok"
fi
. ./.venv/bin/activate

if [ ! -e "./${CENV}" ];then
  echo "conda env not created yet, creating..."
  if [ ! -e ./install_conda.sh ];then
    wget \
      -O install_conda.sh \
      https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh \
    ;
  fi
  bash ./install_conda.sh -b -p "./${CENV}"

  ./${CENV}/bin/conda install -y -n base conda-libmamba-solver
  ./${CENV}/bin/conda create \
    -y \
    --quiet \
    --override-channels \
    --channel conda-forge \
    --channel bioconda \
    --channel defaults \
    --name "${ENV_NAME}" \
    ${MAMBA_SOLVER} \
    sirius-csifingerid=${SIRIUS_VERSION} \
    r-optparse=1.7.1 \
    r-base=4.1.3
  echo "conda env created"
fi

echo ""
echo "===== preparing ====="

oldwd=$(pwd)
tmp=$(mktemp -d)
echo "Working in ${tmp}"
cd "${tmp}"
echo "creating links..."
ln -s "${oldwd}/test-data" test-data
ln -s "${oldwd}/MS2snoop.R" MS2snoop.R
echo "ready to work"

echo ""
echo "===== processing ====="

. "${oldwd}/${CENV}/bin/activate" "${oldwd}/${CENV}/envs/${ENV_NAME}" ;

Rscript ./MS2snoop.R \
  -c ./test-data/compounds_pos.csv \
  -f ./test-data/peaklist_fragments.tsv \
  -p ./test-data/peaklist_precursors.tsv \
  -o ${tmp}/out \
;

echo ""
echo "Error code: ${?}" ;

lines=$(diff ${tmp}/out test-data/compound_fragments_result.txt 2>&1)

echo ""
echo "===== results ====="

if [ "${lines}" = "" ];then
  echo "Result equal to expected."
else
  echo "Some lines are different:"
  echo "${lines}"
fi

echo ""
echo "===== cleaning ====="

echo "Removing ${tmp}..."
rm -rf "${tmp}"

echo ""
echo "Done."
