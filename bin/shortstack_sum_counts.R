#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(readr)
library(tidyr)

option_list <- list(
    make_option(c("-g", "--gff"           ), type="character", default=NULL, metavar="character", help="ShortStack GFF3 file"                        ),
    make_option(c("-c", "--counts"        ), type="character", default=NULL, metavar="character", help="ShortStack counts file"                      ),
    make_option(c("-f", "--count_col"     ), type="character", default=NULL, metavar="character", help="Name of column containing sample count data."),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''  , metavar="character", help="Suffix to remove from sample names."         )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$gff) || is.null(opt$counts)) {
    print_help(opt_parser)
    stop("GFF and counts files must be provided", call.=FALSE)
}

# Read counts
counts           <- read.delim(file=opt$counts, header=TRUE, row.names=NULL, check.names=FALSE, comment.char="#")
colnames(counts) <- gsub(opt$sample_suffix, "", colnames(counts))
colnames(counts) <- gsub(pattern='\\.$', replacement='', colnames(counts))

# Read GFF3
gff <- read_tsv(opt$gff, col_names = FALSE, comment = "#", col_types = cols(.default = "c"))
gff <- gff %>%
    mutate(Locus = sub("^.*?=", "", sub(";.*", "", X9))) %>%
    rename(FeatureType = X3) %>%
    select(Locus, FeatureType)

## Join tables
counts <- merge(gff, counts, by.x = "Locus", by.y = "Name")
counts <- counts[, c("FeatureType", gsub(opt$sample_suffix, "", opt$count_col))]

# Sum counts
res <- counts %>%
    group_by(FeatureType) %>%
    summarise(across(everything(), \(x) sum(as.numeric(x), na.rm = TRUE))) %>%
    rename(biotype = FeatureType)

# Write output
outprefix <- gsub(opt$sample_suffix, "", opt$count_col)
write.table(res, file = paste(outprefix, ".shortstack.counts.vals.txt", sep=""),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
