#!/bin/bash

TISSUE="hippocampus"
EQTL="../../_m/Brainseq_LIBD.allpairs.txt.gz"

## Estimate priors for fine-mapping
/ceph/opt/torus/src/torus.static \
    -d $EQTL --fastqtl -dump_prior priors
