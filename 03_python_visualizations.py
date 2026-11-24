import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import mygene

# 1. Load your cleaned DEG table
df = pd.read_csv(r"D:\Ythdc1_m6A_RNAseq_Portfolio\data\Significant_DEGs_Ythdc1_KO.csv")

# Ensure numeric columns are strictly clean
for col in ["value_1", "value_2", "log2FC", "q_val"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["value_1", "value_2", "log2FC", "q_val"])
df = df[(df["value_1"] >= 0) & (df["value_2"] >= 0)].copy()

# 2. Map RefSeq IDs (NM_...) to official Mouse Gene Symbols
df["refseq_clean"] = df["gene_id"].astype(str).str.replace(r"\..*", "", regex=True)

mg = mygene.MyGeneInfo()
query_res = mg.querymany(
    df["refseq_clean"].unique().tolist(),
    scopes="refseq",
    fields="symbol",
    species="mouse",
    verbose=False
)

symbol_map = {item["query"]: item.get("symbol") for item in query_res if "symbol" in item}
df["Gene_Symbol"] = df["refseq_clean"].map(symbol_map).fillna(df["refseq_clean"])

# Keep unique genes sorted by significance and fold change
annotated = (
    df.sort_values(by=["q_val", "log2FC"], ascending=[True, False])
    .drop_duplicates(subset=["Gene_Symbol"])
    .reset_index(drop=True)
)

# FIGURE 1: Diverging Bar Chart of Top 10 Up & Top 10 Downregulated Genes
top_up = annotated[annotated["log2FC"] > 0].nlargest(10, "log2FC")
top_down = annotated[annotated["log2FC"] < 0].nsmallest(10, "log2FC")
top20 = pd.concat([top_up, top_down]).sort_values("log2FC")

plt.figure(figsize=(9, 6.5), dpi=300)
colors = ["#2b8cbe" if x < 0 else "#e34a33" for x in top20["log2FC"]]

plt.barh(top20["Gene_Symbol"], top20["log2FC"], color=colors, edgecolor="black", alpha=0.85)
plt.axvline(0, color="black", linewidth=1)
plt.xlabel("Log2 Fold Change (Ythdc1 KO / WT)", fontsize=12, fontweight="bold")
plt.ylabel("Mouse Gene Symbol", fontsize=12, fontweight="bold")
plt.title("Top 20 Differentially Expressed Genes in Ythdc1 KO mESCs", fontsize=13, fontweight="bold", pad=12)
plt.grid(axis="x", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(r"D:\Ythdc1_m6A_RNAseq_Portfolio\results\Python_Top20_DEGs_BarChart.png", dpi=300)
plt.close()

# FIGURE 2: Expression Heatmap (Top 12 Upregulated + Top 12 Downregulated Genes)
top24 = pd.concat([
    annotated[annotated["log2FC"] > 0].head(12),
    annotated[annotated["log2FC"] < 0].head(12)
]).reset_index(drop=True)

heatmap_data = pd.DataFrame({
    "WT_mESC": np.log2(top24["value_1"].to_numpy() + 1),
    "Ythdc1_KO_mESC": np.log2(top24["value_2"].to_numpy() + 1)
}, index=top24["Gene_Symbol"].tolist())

plt.figure(figsize=(6.5, 9), dpi=300)
sns.heatmap(
    heatmap_data,
    cmap="RdYlBu_r",
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "Log2(FPKM + 1)"}
)
plt.title("Expression Heatmap: Top 24 Significant Genes (WT vs KO)", fontsize=12, fontweight="bold", pad=12)
plt.ylabel("Mouse Gene Symbol", fontsize=11, fontweight="bold")
plt.tight_layout()
plt.savefig(r"D:\Ythdc1_m6A_RNAseq_Portfolio\results\Python_Top24_Expression_Heatmap.png", dpi=300)
plt.close()

print("Saved both Python plots cleanly! Rows plotted:", len(heatmap_data))
