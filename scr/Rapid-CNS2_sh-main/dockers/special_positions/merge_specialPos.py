#!/usr/bin/env python

import pandas as pd
import glob 

files = glob.glob("/run/special_pos/*.vcf")

names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

mutation={name:name.split("_")[-1][:-4] for name in files}

data = pd.read_csv(files[0], sep="\t", names=names)



data["mutation"]=mutation[files[0]]

for file in files[1:]:
    helper = pd.read_csv(file, sep="\t", names=names)
    helper["mutation"] = mutation[file]
    data = pd.concat([data, helper])
    
data.reset_index(inplace=True)
info_names = list(data.loc[1, "FORMAT"].split(":"))
data[info_names] = data["SAMPLE"].str.split(":", expand=True)

data["F"] = data["ADF"].str.split(",").str[:2].str.join(",")
data[["RF", "AF"]] = data["F"].str.split(",", expand=True)
data["R"] = data["ADR"].str.split(",").str[:2].str.join(",")
data[["RR", "AR"]] = data["R"].str.split(",", expand=True)

data.drop(data.columns.difference(["RF", "AF", "RR", "AR", "DP", "mutation", "CHROM", "POS", "REF", "ALT"]), axis=1, inplace=True)
data = data[["CHROM", "POS", "mutation", "REF", "ALT", "DP", "RF", "AF", "RR", "AR"]]
data.to_excel("/run/special_pos/summary.xlsx", index=False)