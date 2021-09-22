#!/bin/bash

TISSUE="caudate"
PRIOR=(priors/*.prior)

## Running DAP-G
mkdir dapg_out
printf '%s\n' "${PRIOR[@]}" | \
    parallel --jobs 48 "/ceph/opt/dap/dap_src/dap-g -d ${TISSUE}/{/.}.sbams.dat \
                        -p {} -ld_control 0.5 --all -t 2 > dapg_out/{/.}.out"
