for (package in c('optparse', 'rmarkdown','kableExtra','knitr', 'ggplot2')) {
    library(package, character.only=T)
}

#Parse arguments
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix", metavar="character"),
  make_option(c("-i", "--indels"), type="character", default=NULL, 
              help="indel table", metavar="character"),
  make_option(c("-m", "--snv"), type="character", default=NULL, 
              help="snv table", metavar="character"),
  make_option(c("-f", "--fusions"), type="character", default=NULL, 
              help="fusion table", metavar="character"),
  make_option(c("-c", "--cnv_plot"), type="character", default=NULL,
              help="cnv plot", metavar="character"),
  make_option(c("-r", "--rf_details"), type="character", default=NULL,
              help="RF details tsv", metavar="character"),
  make_option(c("-v", "--votes"), type="character", default=NULL,
              help="votes file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-n", "--library"), type="character", default=NULL, 
              help="library", metavar="character"),
  make_option(c("-e", "--coverage"), type="character", default=NULL,
              help="coverage summary", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample", metavar="character"),
  make_option(c("-t", "--mgmt"), type="character", default=NULL,
              help="mgmt prediction", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
prefix <- opt$prefix
indel <- opt$indels
snv <- opt$snv
fusions <- opt$fusions
rf_details <- opt$rf_details
cnv_plot <- opt$cnv_plot
votes <- opt$votes
coverage <- opt$coverage
library <- opt$library
mgmt <- opt$mgmt
sample <- opt$sample


render("rdata/Rapid_CNS2_report_UKHD.Rmd", 
       output_format = "html_document", 
       output_dir=opt$output_dir,
       output_file=paste0(prefix,"_Rapid-CNS2_report.html")
)
