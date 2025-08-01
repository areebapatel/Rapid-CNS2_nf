
# Load required packages with error checking
required_packages <- c('optparse', 'rmarkdown', 'kableExtra', 'knitr')

for (package in required_packages) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Required R package '", package, "' is not installed. Please ensure all required packages are installed in the container."))
  }
}


#Parse arguments
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix", metavar="character"),
  make_option(c("-m", "--mutations"), type="character", default=NULL, 
             help="mutations table", metavar="character"),
  make_option(c("-c", "--cnv_plot"), type="character", default=NULL,
              help="cnv plot", metavar="character"),
  make_option(c("-r", "--rf_details"), type="character", default=NULL,
              help="RF details tsv", metavar="character"),
  make_option(c("-v", "--votes"), type="character", default=NULL,
              help="votes file", metavar="character"),
  make_option(c("-n", "--patient"), type="character", default=NULL, 
              help="patient", metavar="character"),
  make_option(c("-e", "--coverage"), type="character", default=NULL,
              help="coverage summary", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample", metavar="character"),
  make_option(c("-t", "--mgmt"), type="character", default="false", 
             help="mgmt prediction", metavar="character"),
  make_option(c("-a", "--methylartist"), type="character", default="false",
              help="methylartist mgmt plot", metavar="character"),
  make_option(c("-b", "--promoter_mgmt_coverage"), type="integer", default=NULL,
              help="average coverage at mgmt promoter", metavar="character"),
  make_option(c("-g", "--igv_report"), type="character", default="false",
              help="IGV-report html output", metavar="character"),
  make_option(c("-i", "--seq"), type="character", default="Unknown",
             help="Platform used to sequencing; F=MinION/GridION, P=PromethION", metavar="character"),
  make_option(c("-k", "--report_UKHD"), type="character", default=NULL,
              help="R Markdown template file", metavar="character"),
  make_option(c("-l", "--software_ver"), type="character", default=NULL,
              help="Software version", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
prefix <- opt$prefix
mutations <- opt$mutations
rf_details <- opt$rf_details
cnv_plot <- opt$cnv_plot
votes <- opt$votes
coverage <- opt$coverage
patient <- opt$patient
seq <- opt$seq
sample <- opt$sample
mgmt <- opt$mgmt
methylartist_plot <- opt$methylartist
cov <- opt$promoter_mgmt_coverage
igv_report <- opt$igv_report
report_UKHD <- opt$report_UKHD
software_ver <- opt$software_ver

# Check if MGMT file exists
mgmt_status = "false"
if (file.exists(opt$mgmt)) {
    mgmt_status = "true"
}

# Check if methylartist plot exists
methylartist_status = "false"
if (file.exists(opt$methylartist)) {
    methylartist_status = "true"
}

# Define logo files
logo_file <- "rapid_cns2_logo.png"
institution_logo <- "ukhd_logo.png"

# generate the report

inc_igvreport = FALSE
exc_igvreport = TRUE
# lite version
render(report_UKHD, 
       output_format = "html_document", 
       output_file = paste0(prefix,"_Rapid-CNS2_report_lite.html"),
       params = list(logo_file = logo_file, institution_logo = institution_logo))

inc_igvreport = TRUE
exc_igvreport = FALSE
# full version
render(report_UKHD,
       output_format = "html_document",
       output_file = paste0(prefix,"_Rapid-CNS2_report_full.html"),
       params = list(logo_file = logo_file, institution_logo = institution_logo))
