#!/usr/bin/env Rscript
# Writes the AlloppNET BEAST XML with Graham Jones's AlloppDT scripts (not part of this repository).
# Usage: Rscript generate_beast_xml.r <alloppdt_dir> <input_dir> <output_dir> <chain_length> <log_every> <screen_log_every>
# input_dir holds <gene>.nex alignments and taxa_table.txt. alloppnet.XML is written to output_dir and BEAST
# writes its samples (sampledmultrees.txt, sampledparams.txt, sampledgtrees*, DBUGTUNE.txt) there too.

library(ape)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("usage: Rscript generate_beast_xml.r <alloppdt_dir> <input_dir> <output_dir> <chain_length> <log_every> <screen_log_every>")
}
for (f in c("AlloppDT_5beastxml_toplevel.r", "AlloppDT_6beastxml_bits.r")) {
    path <- file.path(args[1], f)
    if (!file.exists(path)) stop(path, " not found, set tools.alloppdt to the folder with the AlloppDT R scripts")
    source(path)
}

# AlloppDT pastes file names straight onto data.dpath, so it needs the trailing slash
input_dir <- paste0(normalizePath(args[2], mustWork = TRUE), "/")
output_dir <- normalizePath(args[3], mustWork = TRUE)
if (!file.exists(file.path(input_dir, "taxa_table.txt"))) stop("taxa_table.txt not found in ", input_dir)

# make.beastxml.AlloppDT.forRealData reads the global alignmentnames
alignmentnames <- sub("\\.nex$", "", list.files(input_dir, pattern = "\\.nex$"))
if (length(alignmentnames) == 0) stop("no .nex files in ", input_dir)

log_every <- args[5]
make.beastxml.AlloppDT.forRealData(list(
    data.dpath = input_dir,
    fpath.taxatable = "taxa_table.txt",
    alignmentnames = alignmentnames,
    beastXMLfilename = "alloppnet.XML",
    beast.chain.length = args[4],
    beast.screen.logevery = args[6],
    beast.params.logevery = log_every,
    beast.gtrees.logevery = log_every,
    beast.multree.logevery = log_every,
    beast.dbugtune.logevery = "10000",
    sampledgtrees.fpathbase = file.path(output_dir, "sampledgtrees"),
    sampledmultrees.fpath = file.path(output_dir, "sampledmultrees.txt"),
    sampledparams.fpath = file.path(output_dir, "sampledparams.txt"),
    DBUGTUNE.fpath = file.path(output_dir, "DBUGTUNE.txt")
))

# AlloppDT writes the XML into data.dpath
if (!file.copy(file.path(input_dir, "alloppnet.XML"), file.path(output_dir, "alloppnet.XML"), overwrite = TRUE)) {
    stop("could not copy alloppnet.XML to ", output_dir)
}
file.remove(file.path(input_dir, "alloppnet.XML"))
