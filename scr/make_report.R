
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
  make_option(c("-b", "--promoter_mgmt_coverage"), type="double", default=NULL,
              help="average coverage at mgmt promoter", metavar="character"),
  make_option(c("-g", "--igv_report"), type="character", default="false",
              help="IGV-report html output", metavar="character"),
  make_option(c("-l", "--report_HTML"), type="character", default="scr/Rapid_CNS2_report_UKHD_HTML.Rmd",
              help="R Markdown template file", metavar="character"),
  make_option(c("-f", "--fusions"), type="character", default="false",
              help="directory of per-caller reportable fusion tables", metavar="character"),
  make_option(c("-y", "--mnpflex"), type="character", default="false",
              help="MNP-Flex predictions table", metavar="character"),
  make_option(c("-w", "--savana_cnv_plot"), type="character", default="false",
              help="SAVANA absolute copy number plot", metavar="character"),
  make_option(c("-u", "--purity_ploidy"), type="character", default="false",
              help="SAVANA purity/ploidy fit", metavar="character"),
  make_option(c("-x", "--egfrviii"), type="character", default="false",
              help="directory of per-caller EGFRvIII checks", metavar="character"),
  make_option(c("-o", "--software_ver"), type="character", default=NULL,
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
sample <- opt$sample
mgmt <- opt$mgmt
methylartist_plot <- opt$methylartist
cov <- opt$promoter_mgmt_coverage
igv_report <- opt$igv_report
software_ver <- opt$software_ver

# --- fusions: merge the per-caller tables, recording which callers agree ------
read_fusions <- function(dir) {
    if (is.null(dir) || dir == "false" || !dir.exists(dir)) return(NULL)
    files <- list.files(dir, pattern = "_fusions_reportable\\.tsv$", full.names = TRUE)
    if (!length(files)) return(NULL)
    rows <- do.call(rbind, lapply(files, function(f) {
        d <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(d) || !nrow(d)) return(NULL)
        d$caller <- sub("^.*_([a-z0-9]+)_fusions_reportable\\.tsv$", "\\1", basename(f))
        d
    }))
    if (is.null(rows) || !nrow(rows)) return(NULL)
    # collapse the same event seen by several callers into one row
    key <- paste(rows$fusion, rows$chrom1, rows$pos1, rows$chrom2, rows$pos2)
    agg <- do.call(rbind, lapply(split(rows, key), function(g) {
        g1 <- g[1, , drop = FALSE]
        g1$caller <- paste(sort(unique(g$caller)), collapse = ", ")
        g1$reads  <- suppressWarnings(max(as.numeric(g$reads), na.rm = TRUE))
        g1
    }))
    agg[order(agg$tier, -as.numeric(agg$reads)), , drop = FALSE]
}
fusions_df <- read_fusions(opt$fusions)

mnpflex <- NULL; mnpflex_qc <- NA_character_; mnpflex_clf <- NA_character_
if (!is.null(opt$mnpflex) && opt$mnpflex != "false" && file.exists(opt$mnpflex)) {
    hdr <- grep("^#", readLines(opt$mnpflex, warn = FALSE), value = TRUE)
    grab <- function(k) { h <- grep(paste0("^# ", k), hdr, value = TRUE)
                          if (length(h)) sub("^[^\t]*\t", "", h[1]) else NA_character_ }
    mnpflex_clf <- grab("classifier"); mnpflex_qc <- grab("qc")
    mnpflex <- tryCatch(read.delim(opt$mnpflex, comment.char = "#",
                                   stringsAsFactors = FALSE), error = function(e) NULL)
}

savana_cnv_plot <- opt$savana_cnv_plot
if (is.null(savana_cnv_plot) || savana_cnv_plot == "false" ||
    !file.exists(savana_cnv_plot)) savana_cnv_plot <- NULL

purity_ploidy <- NULL
if (!is.null(opt$purity_ploidy) && opt$purity_ploidy != "false" &&
    file.exists(opt$purity_ploidy)) {
    purity_ploidy <- tryCatch(read.delim(opt$purity_ploidy, stringsAsFactors = FALSE),
                              error = function(e) NULL)
}

egfrviii_txt <- NA_character_
if (!is.null(opt$egfrviii) && opt$egfrviii != "false" && dir.exists(opt$egfrviii)) {
    ef <- list.files(opt$egfrviii, pattern = "_EGFRvIII\\.txt$", full.names = TRUE)
    hits <- unlist(lapply(ef, function(f) grep("CANDIDATE", readLines(f, warn = FALSE), value = TRUE)))
    egfrviii_txt <- if (length(hits)) paste(hits, collapse = "; ") else "not detected"
}

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



# Set variables for Rmd environment
nextflow_ver <- software_ver

# Set additional variables that Rmd files expect
show_mgmt <- mgmt_status == "true" && methylartist_status == "true"
mgmt_too_low <- !file.exists(mgmt) && !file.exists(methylartist_plot)

# generate the report

inc_igvreport = FALSE
exc_igvreport = TRUE
# lite version - HTML (no embedded IGV report)
render(opt$report_HTML,
       output_format = "html_document",
       output_file = paste0(prefix,"_Rapid-CNS2_report_lite.html"))

# full version - HTML. Only embed the IGV report if it was actually produced.
inc_igvreport = file.exists(igv_report)
exc_igvreport = !inc_igvreport
if (!inc_igvreport) {
    message("No IGV report supplied - the full report will omit the IGV section.")
}
render(opt$report_HTML,
       output_format = "html_document",
       output_file = paste0(prefix,"_Rapid-CNS2_report_full.html"))
