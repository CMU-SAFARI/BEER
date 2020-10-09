#!/usr/bin/env bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd ${SCRIPT_DIR}/z3-4.8.9
python scripts/mk_make.py --prefix=${SCRIPT_DIR}/z3-4.8.9
cd build
make -j 8
make install
