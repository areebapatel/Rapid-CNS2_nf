#!/usr/bin/env Rscript
# Genome-wide absolute copy number from the SAVANA segmentation, with the
# ploidy reference line and CNS driver genes annotated.
suppressPackageStartupMessages({library(ggplot2); library(optparse)})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--segments", type = "character"),
  make_option("--purity_ploidy", type = "character", default = "false"),
  make_option("--genes", type = "character"),
  make_option("--sample", type = "character", default = "sample"),
  make_option("--out", type = "character", default = "cnv.png")
)))

seg <- read.delim(opt$segments, stringsAsFactors = FALSE)
seg <- seg[!is.na(seg$copyNumber) & seg$copyNumber != "", ]
seg$copyNumber <- as.numeric(seg$copyNumber)

ploidy <- NA
if (opt$purity_ploidy != "false" && file.exists(opt$purity_ploidy)) {
  pp <- read.delim(opt$purity_ploidy, stringsAsFactors = FALSE)
  ploidy <- as.numeric(pp$ploidy[1]); purity <- as.numeric(pp$purity[1])
}

chroms <- paste0("chr", c(1:22, "X"))
seg <- seg[seg$chromosome %in% chroms, ]
seg$chromosome <- factor(seg$chromosome, levels = chroms)

# cumulative genome coordinates
len <- tapply(seg$end, seg$chromosome, max)
len[is.na(len)] <- 0
off <- c(0, cumsum(as.numeric(len))[-length(len)])
names(off) <- chroms
seg$gstart <- seg$start + off[as.character(seg$chromosome)]
seg$gend   <- seg$end   + off[as.character(seg$chromosome)]

# call state relative to ploidy
ref <- if (!is.na(ploidy)) ploidy else median(seg$copyNumber, na.rm = TRUE)
seg$state <- ifelse(seg$copyNumber < ref * 0.75, "Loss",
             ifelse(seg$copyNumber > ref * 1.25, "Gain", "Neutral"))

# CNS CNV-relevant genes: the same curated set cnvAnnotated uses for CNVpytor,
# so both copy-number views are annotated against one gene list.
rg <- read.delim(opt$genes, header = FALSE, stringsAsFactors = FALSE)[, 1:4]
names(rg) <- c("chromosome", "start", "end", "gene")
rg <- rg[rg$chromosome %in% chroms, ]
rg$pos  <- (rg$start + rg$end) / 2
rg$gpos <- rg$pos + off[as.character(rg$chromosome)]
rg$cn <- sapply(seq_len(nrow(rg)), function(i) {
  s <- seg[seg$chromosome == rg$chromosome[i] &
           seg$start <= rg$pos[i] & seg$end >= rg$pos[i], ]
  if (nrow(s)) s$copyNumber[1] else NA
})

ymax <- min(10, max(8, ceiling(quantile(seg$copyNumber, 0.995, na.rm = TRUE))))
seg$capped <- seg$copyNumber > ymax
mids <- off + as.numeric(len) / 2

p <- ggplot(seg) +
  geom_rect(aes(xmin = gstart, xmax = gend, ymin = pmin(copyNumber, ymax) - 0.06,
                ymax = pmin(copyNumber, ymax) + 0.06, fill = state)) +
  geom_vline(xintercept = off[-1], colour = "#e6e6e6", linewidth = 0.3) +
  scale_fill_manual(values = c(Loss = "#0047b9", Neutral = "#b0b8c1", Gain = "#ff5c00")) +
  scale_x_continuous(breaks = mids, labels = sub("chr", "", chroms), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, ymax + 5), breaks = 0:ymax) +
  labs(x = NULL, y = paste0("absolute copy number (capped at ", ymax, ", \u25b3 = above)"),
       title = paste0(opt$sample,
                      if (!is.na(ploidy)) sprintf("   (purity %.0f%%, ploidy %.2f)",
                                                  purity * 100, ploidy) else "")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.line = element_line(colour = "#ccd6e5"),
        axis.ticks = element_line(colour = "#ccd6e5"),
        panel.grid.major.y = element_line(colour = "#f0f2f5"),
        plot.title = element_text(size = 11, colour = "#0047b9"))

if (!is.na(ploidy))
  p <- p + geom_hline(yintercept = ploidy, linetype = "dashed", colour = "#28348a", linewidth = 0.4)

rg <- rg[!is.na(rg$cn), ]
rg <- rg[order(rg$gpos), ]
if (nrow(rg)) {
  rg$capped <- rg$cn > ymax
  rg$y      <- pmin(rg$cn, ymax)
  # alternate the label offset when two genes are close on the x axis
  gap  <- diff(c(-Inf, rg$gpos))
  near <- gap < (max(seg$gend) * 0.02)
  lev  <- 0; offs <- numeric(nrow(rg))
  for (i in seq_len(nrow(rg))) { lev <- if (near[i]) lev + 1 else 0; offs[i] <- lev %% 3 }
  rg$ylab <- rg$y + 0.45 + offs * 1.6
}
if (nrow(rg))
  p <- p +
    geom_point(data = rg[!rg$capped, ], aes(x = gpos, y = y),
               size = 1.4, colour = "#c61826") +
    # open upward triangle: the true copy number is above the axis limit
    geom_point(data = rg[rg$capped, ], aes(x = gpos, y = y),
               size = 2.2, shape = 2, stroke = 0.8, colour = "#c61826") +
    geom_text(data = rg, aes(x = gpos, y = ylab, label = gene),
              size = 2.6, angle = 90, hjust = 0, colour = "#c61826")

ggsave(opt$out, p, width = 13, height = 5, dpi = 150)
cat("wrote", opt$out, "with", nrow(seg), "segments and", nrow(rg), "gene labels\n")
