# Methylation Classification
- ad-hoc random forest classifier, trained on 450k methylation array reference data set, HD classifier v11
- top overlapping 100k sites are chosen (mean decrease in accuracy)
- variance filter for selecting top 10k probes
- random forest with 20k trees is trained, recalibration by one vs all GLM for each class 
- ground truth is inferred from EPIC array data (FFPE samples) with Hd classifier predictions v11b4
- families are inferred aggregating over methylation subclasses from reference set

## Docker Versions
- v1.0 initial set-up as coded by Areeba