#!/usr/bin/env python

import sys
import pandas as pd 
import re
PREFIX=sys.argv[1]

try: 
	deepvar = pd.read_csv(sys.argv[2])
except FileNotFoundError:
	print("no deepvariant vcf")


if len(sys.argv) == 1:
	sys.exit("no variants files")

df = deepvar

df.drop("cytoBand", axis=1, inplace=True)

df.rename({"Chr":"chr", "Start":"start", "End":"end", "Ref":"ref", "Alt":"alt", "Func.refGene":"type", "Gene.refGene":"gene", "AAChange.refGene":"transcript", "1000g2015aug_all":"1000G", "cosmic68":"COSMIC", "ExonicFunc.refGene":"effect"}, inplace=True, axis=1)

df["1000G"] =df["1000G"].fillna(0)
df = df[df["1000G"] <= 0.01]
df["1000G"] = df["1000G"].apply(lambda x: round(x, 2))

# for now
helper = {"Chr":"chr", "Start":"start", "End":"end", "Ref":"ref", "Alt":"alt", "Func.refGene":"type", "Gene.refGene":"gene", "AAChange.refGene":"transcript", "X1000g2015aug_all":"1000G", "ExonicFunc.refGene":"effect"}

df = df[list(helper.values())+["VAF", "DP"]]
df = df.loc[:,~df.columns.duplicated()].copy()
try:
    df["transcript"] = df["transcript"].str.split(",").str[0]
except Exception:
    pass
df["transcript"] = df["transcript"].replace(r"exon[0-9]{1,2}:", "", regex=True)

df["dups"] = df.duplicated(['chr', 'start', 'end', 'ref', 'alt'], keep=False)

df.drop("dups", inplace=True, axis=1)

df.drop_duplicates(['chr', 'start', 'end', 'ref', 'alt'], inplace=True)

print(df)
df = df[df["chr"].isin(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr11", "chr10", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"])]

df["pos"] = df["chr"].astype(str)+":"+df["start"].astype(str)

df = df[["pos", "ref", "alt", "type", "gene", "VAF", "DP", "transcript"]] # "1000G", "effect"

df["gene"] = df["gene"].str.split(";").str[0]

df = df[~df["type"].isin(["ncRNA_intronic", "ncRNA_exonic", "upstream", "UTR3"])]

indels = df[~((df["ref"].str.len()==df["alt"].str.len())) | ((df["ref"]=="-") | (df["alt"]=="-"))]
snvs = df[((df["ref"].str.len()==df["alt"].str.len()) & (df["ref"]!="-") & (df["alt"]!="-"))]

indels.to_csv(f"{PREFIX}.indels.csv", index=False)
snvs.to_csv(f"{PREFIX}.snvs.csv", index=False)