"""
Combine parquet files into one. This uses with to reduce memory
by only needing to load one files at a time.
"""
import argparse
import session_info
import pandas as pd
from pathlib import Path

def combine_files(output, label):
    data_dir = Path(output+"/")
    with open("%s.sex_interaction.txt" % label, "a+") as handle:
        for ii, tmpfile in enumerate(data_dir.glob(label+'*')):
            write_header = ii == 0
            pd.read_csv(tmpfile, sep='\t')\
              .to_csv(handle, header=write_header, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=str, default="output")
    parser.add_argument('--label', type=str, default="lfsr")
    args=parser.parse_args()

    # Combine files
    combine_files(args.output, args.label)
    
    # Session Information
    session_info.show()


if __name__ == "__main__":
    main()
