# Do NOT auto-install: that writes into the user's home library at run time,
# which is both a surprise side effect and unreproducible. The container
# provides these packages; if they are missing, fail loudly instead.
for (package in c('optparse', 'reshape2')) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package '", package, "' is not installed. ",
         "Please run this inside the pipeline container.")
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input_file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output_file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

var_file <- read.csv(opt$input)

# Optional ANNOVAR databases: cosmic70 needs a COSMIC licence and CLNSIG comes
# from clinvar. Fill in "." so the filters below behave as "no annotation"
# rather than erroring when those protocols were not run.
for (col in c("cosmic70", "CLNSIG", "avsnp151")) {
  if (!col %in% colnames(var_file)) {
    message("Note: '", col, "' not present in ANNOVAR output - filling with '.'")
    var_file[[col]] <- "."
  }
}

# 1000G EUR allele frequency arrives as character, with "." meaning the variant
# is ABSENT from 1000G - i.e. not a common polymorphism, and exactly what a
# somatic driver looks like. Comparing "." numerically yields NA, and which()
# then drops the row, so novel variants were being silently discarded. Coerce
# once and treat absence as frequency 0 so those variants are retained.
af_eur <- suppressWarnings(as.numeric(var_file$X1000g2015aug_eur))
n_novel <- sum(is.na(af_eur))
af_eur[is.na(af_eur)] <- 0
var_file$af_eur <- af_eur
message("Variants absent from 1000G (treated as novel, AF=0): ", n_novel)

is_rare      <- var_file$af_eur < 0.001
in_cosmic    <- var_file$cosmic70 != "."
is_upstream  <- var_file$Func.refGene == "upstream"
not_synon    <- var_file$ExonicFunc.refGene != "synonymous SNV"

keep <- (is_rare | in_cosmic | is_upstream) & not_synon
no_syn <- var_file[which(keep), ]
no_syn <- no_syn[order(no_syn$af_eur), ]
message("Variants retained after filtering: ", nrow(no_syn), " of ", nrow(var_file))

table <- no_syn[,c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","cytoBand","avsnp151","X1000g2015aug_eur","cosmic70","CLNSIG" )]

# Depth and allele fraction come from the Clair3 genotype column. ANNOVAR keeps
# the whole tab-delimited VCF record in a single Otherinfo field (convert2annovar
# was run with -includeinfo), so split on tabs and locate FORMAT/SAMPLE by name -
# Clair3 emits GT:GQ:DP:AD:AF but the field order is not guaranteed.
gt_field <- function(blob, key) {
    if (is.na(blob)) return(NA_character_)
    f <- strsplit(blob, "\t", fixed = TRUE)[[1]]
    i <- grep("^GT:", f)
    if (!length(i) || i[1] + 1 > length(f)) return(NA_character_)
    keys <- strsplit(f[i[1]],     ":", fixed = TRUE)[[1]]
    vals <- strsplit(f[i[1] + 1], ":", fixed = TRUE)[[1]]
    k <- match(key, keys)
    if (is.na(k) || k > length(vals)) NA_character_ else vals[k]
}
oi_cols <- grep("^Otherinfo", names(no_syn), value = TRUE)
blob <- if (length(oi_cols)) {
    apply(no_syn[, oi_cols, drop = FALSE], 1,
          function(r) { g <- r[grepl("GT:", r, fixed = TRUE)]; if (length(g)) g[1] else NA_character_ })
} else rep(NA_character_, nrow(no_syn))

depth <- vapply(blob, gt_field, character(1), key = "DP", USE.NAMES = FALSE)
vaf   <- vapply(blob, gt_field, character(1), key = "AF", USE.NAMES = FALSE)
vaf   <- ifelse(is.na(vaf), NA_character_,
                paste0(round(suppressWarnings(as.numeric(vaf)) * 100, 1), "%"))

var_list <- strsplit(table$AAChange.refGene[1],split = ":")
alteration <- lapply(table$AAChange.refGene, function(x){
  var_list <- strsplit(x,split = ":")
  alter <- paste(var_list[[1]][3],var_list[[1]][4],strsplit(var_list[[1]][5],split = ",")[[1]][1],sep = ";")
  return(alter)
})

cosmic <- lapply(table$cosmic70, function(x){
  var_list <- strsplit(x,split = ";")
  alter <- strsplit(var_list[[1]][1],split = ",")[[1]][1]
  return(alter)
})

filtered <- cbind("Chr"=table$Chr,"Start"=table$Start,"End"=table$End,"Func"=table$Func.refGene,"Gene"=table$Gene.refGene,"ExonicFunc"=table$ExonicFunc.refGene,"AAChange"=reshape2::melt(alteration)["value"],"cytoBand"=table$cytoBand,"1000g_EUR"=table$X1000g2015aug_eur,"COSMIC"=reshape2::melt(cosmic)["value"], "CLNSIG" = table$CLNSIG, "Coverage"=depth, "VAF"=vaf)

colnames(filtered) <- c("Chr","Start","End","Func","Gene","ExonicFunc","AAChange","cytoBand","1000g_EUR","COSMIC","CLNSIG","Coverage","VAF")


write.csv(filtered,file = opt$output,row.names = F)
