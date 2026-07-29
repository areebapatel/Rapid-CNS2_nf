# Do NOT auto-install: that writes into the user's home library at run time,
# which is both a surprise side effect and unreproducible. The container
# provides these packages; if they are missing, fail loudly instead.
for (package in c('optparse')) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package '", package, "' is not installed. ",
         "Please run this inside the pipeline container.")
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="mgmt bed file", metavar="character"),
  make_option(c("-p", "--probes"), type="character", default=NULL, 
              help="top probes", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="mgmt prediction model", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

load(opt$probes)
load(opt$model)
mgmt_meth <- read.delim(opt$input,header = F)
mgmt_meth <- subset(mgmt_meth, V2 %in% pred_pos)

if (nrow(mgmt_meth) == 0) {
  stop("No MGMT promoter CpG sites from the prediction model were found in ", opt$input)
}

# modkit bedMethyl is all-tab delimited: column 11 is the percent modified
# (column 10 is Nvalid_cov). Older mixed-delimiter output packed both into V10.
mgmt_average <- mean(as.numeric(mgmt_meth$V11), na.rm = TRUE)
mgmt <- data.frame(average=mgmt_average)
pred <- predict(log.model,newdata = mgmt,type = "response")
mgmt$pred <- pred[[1]]
if (pred[[1]] <0.5){mgmt$status <- "Unmethylated"} else {mgmt$status <- "Methylated"}
write.csv(mgmt,file=paste0(opt$out_dir,"/",opt$sample,"_mgmt_status.csv"),row.names = F)