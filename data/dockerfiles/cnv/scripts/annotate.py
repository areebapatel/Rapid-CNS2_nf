#!/usr/bin/env python

import sys
import pandas as pd 

cnv = pd.read_csv(sys.argv[1], sep="\t", names=["type", "pos", "size", "cnv", "p_val", "p_val_2", "p_val_3", "p_val_4", "Q0", "pN", "dG"])

cnv["chr"] = cnv["pos"].str.split(":").str[0]
cnv["start"] = cnv["pos"].str.split(":").str[1].str.split("-").str[0].astype(int)
cnv["end"] = cnv["pos"].str.split(":").str[1].str.split("-").str[1].astype(int)
bed = pd.read_csv("/genes.bed", sep="\t", names=["chr", "start", "end", "annotation"])

chr_of_interest = set(bed["chr"].tolist())

cnv = cnv[cnv["chr"].isin(chr_of_interest)]

typis = {0:"full", 1:"partial"}

alterations=pd.DataFrame(columns=["gene", "genomic location", "alteration", "type", "position", "size"])
index_alt=0
for index, (chr, start, end, anno) in bed.iterrows():
    tmp = cnv[cnv["chr"] == chr]
    
    # gene complete in bin
    tmp1 = tmp[(tmp["start"] <= start) & (end <= tmp["end"])]
    tmp1["coverage"]=0
    print(anno, start, end)
    
    # start of gene in bin
    tmp2 = tmp[(tmp["start"] <= start) & (start <= tmp["end"])]
    tmp2["coverage"]=1
    
    # end of gene in bin
    tmp3 = tmp[(tmp["start"] <= end) & (end <= tmp["end"])]
    tmp3["coverage"]=1
    
    covered = pd.concat([tmp1, tmp2, tmp3])
    covered = covered.sort_values("coverage")
    covered = covered.drop_duplicates(["type", "pos"])
    covered = covered.reset_index()
    
    if len(covered) == 0:
        pass
    else:
        for index in covered.index.tolist():

            alterations.at[index_alt, "gene"] = anno
            alterations.at[index_alt, "genomic location"] = f"{chr}:{start}-{end}"
            alterations.at[index_alt, "alteration"] = covered.at[index, "type"]
            alterations.at[index_alt, "position"] = covered.at[index, "pos"]
            alterations.at[index_alt, "size"] = int(covered.at[index, "size"])
            alterations.at[index_alt, "type"] = typis[covered.at[index, "coverage"]]
            index_alt+=1

alterations = alterations.sort_values("gene")
alterations.to_excel(sys.argv[2], index=False)