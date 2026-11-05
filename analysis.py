import pandas as pd # pandas 3.x, окружение bioinf
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

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



groups = {
    "GSM4232567": "control", "GSM4232568": "control", "GSM4232569": "control",
    "GSM4232570": "NaIO3", "GSM4232571": "NaIO3", "GSM4232572": "NaIO3",
}
colors = {"control": "tab:blue", "NaIO3": "tab:red"}
os.makedirs("figures", exist_ok=True)

samples = list(log_cpm.columns)
fig, ax = plt.subplots(figsize=(8, 5))
bp = ax.boxplot([log_cpm[s] for s in samples], patch_artist=True)
for patch, s in zip(bp["boxes"], samples):
    patch.set_facecolor(colors[groups[s]])
ax.set_xticks(range(1, len(samples) + 1))
ax.set_xticklabels([f"{groups[s]}\n{s[-3:]}" for s in samples])
ax.set_ylabel("log2(CPM + 1)")
ax.set_title("Распределение экспрессии по образцам")
fig.savefig("figures/boxplot.png", dpi=150, bbox_inches="tight")

X = log_cpm.T
pca = PCA(n_components=2)
coords = pca.fit_transform(X)
var = pca.explained_variance_ratio_ * 100
print("Доля разброса PC1 и PC2, %:", var.round(1))

fig, ax = plt.subplots(figsize=(6,5))
for s, (x,y) in zip(X.index, coords):
    ax.scatter(x, y, color=colors[groups[s]], s=80)
    ax.annotate(s[-3:], (x, y), textcoords="offset points", xytext=(6, 6))
ax.set_xlabel(f"PC1 ({var[0]:.0f}%)")
ax.set_ylabel(f"PC2 ({var[1]:.0f}%)")
ax.set_title("PCA: control (синий) и NaIO3 (красный)")
fig.savefig("figures/pca.png", dpi=150, bbox_inches="tight")
