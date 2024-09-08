#!/usr/bin/env python

import sys
import pandas as pd 
import plotly.express as px
import plotly

summary = pd.read_csv(sys.argv[1], sep="\t")
sample =sys.argv[1].split("/")[-1].split(".")[0]

total = summary[~summary["chrom"].str.endswith("_region")]
total["reads"]="total"
region= summary[summary["chrom"].str.endswith("_region")]
region["chrom"] = region["chrom"].str.split("_").str[0]
region["reads"]="on target"

df = pd.concat([total, region])

df.rename({"chrom": "chromosome", "mean": "coverage"}, inplace=True, axis=1)
print(df)



plt = px.scatter(df, x="chromosome", y="coverage", color="reads", color_discrete_sequence=["mediumpurple", "mediumseagreen"], title=f"{sample}: on-target coverage")

plt.update_layout(yaxis_range=[0, max(list(map(int, df["coverage"])))+10], template="none")

plt.update_traces(marker=dict(size=13,
                              line=dict(width=1,
                                        color='White'), 
                              ),
                  selector=dict(mode='markers'))

plt.update_xaxes(gridcolor="#d6d4d4")
plt.update_yaxes(gridcolor="#d6d4d4")


plotly.io.write_image(plt, f"/run/coverage/{sample}.coverage.png", format="png")
