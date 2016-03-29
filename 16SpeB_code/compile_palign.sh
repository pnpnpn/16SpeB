#!/bin/bash

set -x
set -euo pipefail
IFS=$'\n\t'

cd palign/code
make clean
make
cd ../../
