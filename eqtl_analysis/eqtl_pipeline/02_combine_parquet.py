"""
Combine parquet files into one. This uses with to reduce memory
by only needing to load one files at a time.
"""
import pandas as pd
from pathlib import Path

def main():
    data_dir = Path(".")
    with open("BrainSEQ_TOPMed.interaction.txt", "a+") as tgz_handle:
        for i, parquet_file in enumerate(data_dir.glob('*.parquet')):
            write_header = i == 0
            pd.read_parquet(parquet_file)\
              .to_csv(tgz_handle, header=write_header,
                      index=False, sep='\t')


if __name__ == "__main__":
    main()
