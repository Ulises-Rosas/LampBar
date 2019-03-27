#!/usr/bin/Rscript

suppressMessages({
  
  library(seqinr)
  library(optparse)
})

option_list = list(
  
  make_option( 
    c("-a", "--alignment_file"),
    type="character",
    default=NULL,
    help="Alignment"
    )
  )

opt_parser = OptionParser(usage = "%prog [-a ##] ",option_list = option_list, description = "")

opt = parse_args(opt_parser)
# alignment_file = "cificusOnlyAlnNoHyphenAln"

writeLines(paste0("> ", opt$`alignment_file`))
writeLines(
  seqinr::c2s(
    chars = seqinr::consensus(
      method = "majority",
      seqinr::read.alignment(
        opt$`alignment_file`, format = "fasta"
        )
      )
    )
  )

