library(optparse)

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

# filtering sites that are preselected
mgmt_meth <- read.delim(opt$input,header = F)
mgmt_meth <- subset(mgmt_meth, V2 %in% pred_pos)

# average of v11 predicitons
mgmt_average <- mean(mgmt_meth$V11)
mgmt <- data.frame(average=mgmt_average)

pred <- predict(log.model,newdata = mgmt,type = "response")
mgmt$pred <- pred[[1]]

if (pred[[1]] <0.5){mgmt$status <- "Unmethylated"} else {mgmt$status <- "Methylated"}
write.csv(mgmt,file=paste0(opt$out_dir,"/",opt$sample,"_mgmt_status.csv"),row.names = F)