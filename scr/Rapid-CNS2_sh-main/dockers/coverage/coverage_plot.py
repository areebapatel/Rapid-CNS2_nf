#!/usr/bin/env python

import sys
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

summary = pd.read_csv(sys.argv[1], sep="\t")
sample  = sys.argv[1].split("/")[-1].split(".")[0] #getting sample name from input file

total = summary[~summary["chrom"].str.endswith("_region")]
region= summary[summary["chrom"].str.endswith("_region")]
region["chrom"] = region["chrom"].str.split("_").str[0]

sns.set_theme()

COLOR="gray"
mpl.rcParams['text.color'] = COLOR
mpl.rcParams['axes.labelcolor'] = COLOR
mpl.rcParams['xtick.color'] = COLOR
mpl.rcParams['ytick.color'] = COLOR

# plotting on target coverage as well as total coverage
sns.scatterplot(data=total, x="chrom", y="mean", label="total", zorder=3, s=50)
sns.scatterplot(data=region, x="chrom", y="mean", label="on target", zorder=3, s=50)


# defining plot parameters
plt.ylim(0, max(max(list(map(int, total["mean"]))), max(list(map(int, region["mean"]))))+10)
plt.xticks(rotation=90)
plt.legend(title="reads")
plt.xlabel("chromosomes")
plt.ylabel("mean coverage")
plt.grid(zorder=0, axis="y", linestyle="--")

plt.tight_layout()
plt.title(f"{sample}: on-target coverage")

plt.savefig(f"/run/coverage/{sample}.coverage.pdf")
plt.savefig(f"/run/coverage/{sample}.coverage.png")
