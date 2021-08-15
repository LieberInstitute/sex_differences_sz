#!/bin/bash

for FEATURE in genes transcripts exons junctions; do
    Rscript generate_covs.R --feature $FEATURE
done
