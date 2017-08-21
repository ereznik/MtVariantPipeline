#!/usr/bin/env Rscript

'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...) {
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

write_batch_script <- function(samples, out, dir, genome, squish=T) {
  sink(out)
  ids <- unique(samples$patient)
  for (id in ids) {
    cat('new')
    # cat(paste0('\ngenome ', genome))
    dat <- samples[patient == id]
    bams_to_load <- paste(dat$tumor,
                          dat$normal, 
                          sep = ",")
    cat(paste0('\nload ', bams_to_load))
    cat(paste0('\nsnapshotDirectory ', dir))
    cat('\nmaxPanelHeight 600\n')
    for (i in 1:nrow(dat)) {
      cat(paste0('\ngoto chr', dat$chr[i], ':', dat$start[i]))
      cat(paste0('\nsort base'))
      cat(paste0('\ngoto chr', dat$chr[i], ':', dat$start[i]-1,'-',dat$end[i]+1))
      if (squish) {
        cat('\nsquish')
      }
      cat(paste0('\nsnapshot ', dat$patient[i], "_", dat$gene[i], "_", dat$chr[i], "_", dat$start[i], '.png'))
      cat("\n")
    }
  }
  sink()
}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse', 'stringr')
  tmp <- lapply(pkgs, require, character.only = T)
  rm(tmp)
  
  parser=ArgumentParser()
  parser$add_argument('-s', '--samples', type='character', help='Tab-delimited text file mapping sample IDs to bam files')
  parser$add_argument('-o', '--out', type='character', help='Output batch script file')
  parser$add_argument('-d', '--dir', type='character', help='Screenshot directory')
  parser$add_argument('-g', '--genome', type='character', default='b37', help='Genome version')
  parser$add_argument('-sq', '--squish', action="store_true", default=T, help='Squish reads')
  args=parser$parse_args()
  
  samples <- fread(args$samples)
  out <- args$out
  dir <- args$dir
  genome <- args$genome
  #squish <- args$squish # set squish to false for now
  squish <- FALSE
  
  setnames(samples, c('patient', 
                      'gene', 
                      'chr', 
                      'start',
		      'end', 
                      'tumor', 
                      'normal'))
  write_batch_script(samples, out, dir, 
                     genome, squish)
  
}
