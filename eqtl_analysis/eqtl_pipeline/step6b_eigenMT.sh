#!/bin/bash

## Clean intermediate files
rm *gz *txt
## Combine and preform FDR correction
python ../_h/combine_N_fdrcorrection.py
## Remove old files
rm *tsv
