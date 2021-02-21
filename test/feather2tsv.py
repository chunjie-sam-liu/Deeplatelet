import feather
import pandas as pd
import os


def feather2tsv(name):
    filepath = f"data/rda/riskscore/pfs.{name}.test.feather"
    df = feather.read_dataframe(source=filepath)
    df.to_csv(f"data/rda/riskscore/pfs.{name}.test.tsv", index=False, sep="\t")


if __name__ == "__main__":
    feather2tsv("train")
    feather2tsv("val")
    feather2tsv("test1")
    feather2tsv("test2")
