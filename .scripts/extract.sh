#!/usr/bin/env bash

WORKDIR=${INIT_CWD:-$(pwd)}
DIR=${1}

cd $WORKDIR

cp $DIR/{*.xyz,*.inp,*.sh,*.log} .