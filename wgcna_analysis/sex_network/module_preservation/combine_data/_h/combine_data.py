#### Combine data

import session_info
import polars as pl
import pandas as pd

def load_data(tissue, dx):
    new_dx = dx.replace("SZD", "SZ")
    fn = f"../../_m/{tissue.lower()}/{dx}/zsummary_table.tsv"
    return pl.from_pandas(pd.read_csv(fn, sep='\t')\
                          .reset_index()\
                          .rename(columns={"index":"module"}))\
             .with_columns([
                 pl.col("Zsummary.pres").cast(pl.Float64),
                 pl.lit(tissue).alias("region"),
                 pl.lit(new_dx).alias("diagnosis")
             ])


def merge_data():
    dt = pl.DataFrame()
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        for dx in ["CTL", "SZD"]:
            dt = pl.concat([dt, load_data(tissue, dx)])
    return dt


def main():
    ## Generate data
    merge_data().write_csv("module_preservation.byDx.csv")
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
