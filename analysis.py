import pandas as pd # pandas 3.x, оружение bioinf
import numpy as np 

counts_file = "data/GSE142591_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annot_file = "data/Human.GRCh38.p13.annot.tsv.gz"

counts = pd.read_csv(counts_file, sep="\t", index_col=0)
annot = pd.read_csv(annot_file, sep="\t", usecols=["GeneID", "Symbol"])

df = counts.reset_index().merge(annot, on="GeneID")
print("Генов в счетах:", len(counts), "После merge:", len(df))
print(df[df["Symbol"].isin(["HMOX1", "SLC7A11"])])



lib_sizes = counts.sum()
print(lib_sizes)

cpm = counts / lib_sizes * 1_000_000
log_cpm = np.log2(cpm + 1) 

keep = (cpm > 1).sum(axis=1) >= 3
log_cpm = log_cpm[keep]
print("Генов до фильтра:", len(counts), "После:", len(log_cpm))
print(log_cpm.loc[[3162,23657]].round(1))

