#!/bin/sh

set -ex

SRC_DIR=$1
INSTALL_DIR=$2

if [ ! -e ${INSTALL_DIR}/.install_finished ]
then
    echo 'Neuron was not fully installed in previous build, installing ...'
    mkdir -p ${SRC_DIR}
    cd ${SRC_DIR}
    echo "Downloading NEURON ..."
    rm -rf nrn
    git clone --depth 1 https://github.com/nrnhines/nrn.git
    cd nrn
    echo "Preparing NEURON ..."
    ./build.sh
    echo "Configuring NEURON ..."
    ./configure --prefix=${INSTALL_DIR} --without-x --with-nrnpython=python --disable-rx3d
    echo "Building NEURON ..."
    make -j4 >make.log 2>&1
    echo "Installing NEURON ..."
    make -j4 install
    cd src/nrnpython
    python setup.py install
    echo "Testing NEURON import ...."
    python -c 'import neuron'

    touch -f ${INSTALL_DIR}/.install_finished
    echo "NEURON successfully installed"
else
    echo 'Neuron was successfully installed in previous build, not rebuilding'
fi
