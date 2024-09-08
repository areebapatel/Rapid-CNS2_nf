for (package in c('optparse', 'rmarkdown','kableExtra','knitr', 'ggplot2', 'openxlsx')) {
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
              help="mgmt prediction", metavar="character"),
  make_option(c("--tissue"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--seqdate"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--prepkit"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--panel"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--refgen"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--device"), type="character", default=NULL, help="tissue type", metavar="character"),
  make_option(c("--specpos"), type="character", default=NULL, help="special positions", metavar="character"),
  make_option(c("--classplot"), type="character", default=NULL, help="classification plot", metavar="character")

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
tissue <- opt$tissue
seqdate <- opt$seqdate
prepkit <- opt$prepkit
panel <- opt$panel
device <- opt$device
specpos <- opt$specpos
classplot <- opt$classplot

if (panel == "NPHD2022A") {
  num_genes <- 169
} else {
  num_genes <- "NaN"
}

refgen <- opt$refgen

if (refgen == "hg19.fa") {
  refgen <- "HG19"
} else if (refgen == "hg38.fa") {
  refgen <- "HG38"
}


render("rdata/Rapid_CNS2_report_UKHD_pdf.Rmd", 
       output_format = "pdf_document", 
       output_dir=opt$output_dir,
       output_file=paste0(prefix,"_Rapid-CNS2_report.pdf")
)

