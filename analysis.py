import pandas as pd # pandas 3.x, оружение bioinf

counts_file = "data/GSE142591_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annot_file = "data/Human.GRCh38.p13.annot.tsv.gz"

counts = pd.read_csv(counts_file, sep="\t", index_col=0)
annot = pd.read_csv(annot_file, sep="\t", usecols=["GeneID", "Symbol"])

df = counts.reset_index().merge(annot, on="GeneID")
print("Генов в счетах:", len(counts), "После merge:", len(df))
print(df[df["Symbol"].isin(["HMOX1", "SLC7A11"])])

