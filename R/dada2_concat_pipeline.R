#' Unified pipeline for analyzing data with different V-regions in dada2
#'
#' Perform trimming&filtering, denoising, merging, chimera removal and taxonomic assignment steps of dada2.
#' If regions include amplicons larger than trimf+trimr, the justConcatenate option should be set.
#' Otherwise, most reads will not be merged resulting in poor results.
#' @param fastq_path Input folder with fastq read data
#' @param samplefile_path Input file with sample names and corresponding R1 and R2 read files
#' @param out_prefix Prefix for output files
#' @param filter_path Output path for filtered and trimmed fastq files (default: fastq_filtered/)
#' @param silva_path Path to SILVA database
#' @param trimf Fixed trimming length for R1 (default: 150)
#' @param trimr Fixed trimming length for R2 (default: 150)
#' @param trimleft trim at the 5' end of reads (e.g., to remove primers, default: 0)
#' @param justConcatenate Don't merge denoised R1 and R2 reads but concatenate them. (Default: TRUE)
#' @param maxEE max. expected error rate (per 100 basepairs of trimmed sequence) used for filtering (default: 1)
#' @param onlyNfirst only for testing purpose. If set to a number n, only the first n samples will be included in the analysis.
#' @export
#'

dada2_concat_pipeline <- function(fastq_path,
                                  samplefile_path,
                                  out_prefix,
                                  filter_path = "fastq_filtered/",
                                  silva_path = "/data/db/SILVA/dada2/silva_nr_v138_train_set.fa.gz",
                                  trimf = 150,
                                  trimr = 150,
                                  trimleft = 0,
                                  justConcatenate = T,
                                  maxEE = 1,
                                  onlyNfirst = NA)
{

  # required packages
  require(dada2)
  require(phyloseq)
  require(reshape2)

  # define helper function
  getN <- function(x) sum(getUniques(x))

  # further parameters
  seqtab_out <- paste0(out_prefix, ".seqtab.rds") # Output nochim sequence table
  taxtab_out <- paste0(out_prefix, ".taxtab.rds") # Output taxonomy table
  track_out  <- paste0(out_prefix, ".track.rds")  # track reads

  maxEE <- maxEE * (trimf+trimr)/200  # max. expected error

  # add trimleft to trim length
  if(length(trimleft) == 1){
    trimf = trimf + trimleft
    trimr = trimr + trimleft
  } else {
    trimf = trimf + trimleft[1]
    trimr = trimr + trimleft[2]
  }

  errorbase = 1e8 # number of bases used for error learning


  # load required data
  if(!is.na(onlyNfirst)){
    filetable <- read.table(samplefile_path)[1:onlyNfirst,]
  } else {
    filetable <- read.table(samplefile_path)
  }
  colnames(filetable) <- c("sample_id", "F_reads", "R_reads")
  F_raw     <- paste0(fastq_path, as.character(filetable[,2]))
  R_raw     <- paste0(fastq_path, as.character(filetable[,3]))
  sample_id <- as.character(filetable[,1])

  # set filter output file names
  F_filtered <- paste0(filter_path, sample_id, "_F.fastq.gz")
  R_filtered <- paste0(filter_path, sample_id, "_R.fastq.gz")
  names(F_filtered) <- sample_id
  names(R_filtered) <- sample_id

  # filter and trim
  filtres <- filterAndTrim(F_raw, F_filtered,
                           R_raw, R_filtered,
                           truncLen = c(trimf,trimr),
                           trimLeft = trimleft,
                           maxN = 0,
                           maxEE = c(maxEE,maxEE),
                           rm.phix = T,
                           compress = T,
                           multithread = T)
  # TODO: export filtres

  # prepare derep objects
  F_derep <- derepFastq(F_filtered, n = 2e5)
  R_derep <- derepFastq(R_filtered, n = 2e5)

  # learn error models
  F_err <- learnErrors(F_derep, nbases = errorbase, multithread = T, randomize = T, MAX_CONSIST = 20)
  R_err <- learnErrors(R_derep, nbases = errorbase, multithread = T, randomize = T, MAX_CONSIST = 20)

  # sample inference, using pseudo-pooling
  F_dada <- dada(F_derep, err = F_err, pool = "pseudo", multithread = T)
  R_dada <- dada(R_derep, err = R_err, pool = "pseudo", multithread = T)

  # concatenate F and R reads
  merged <- mergePairs(F_dada, F_derep, R_dada, R_derep, justConcatenate = justConcatenate)

  seqtab <- makeSequenceTable(merged)

  # remove chimeras
  nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=T, verbose = T)

  # assign taxonomy
  tax <- assignTaxonomy(nochim, refFasta = silva_path, multithread = T)
  colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  # track number of reads throughout process
  if(nrow(filtres)>1){
    trackReads <- cbind(filtres,
                   sapply(F_dada, getN),
                   sapply(R_dada, getN),
                   sapply(merged, getN),
                   rowSums(nochim))
  } else {
    trackReads <- cbind(filtres,
                   getN(F_dada),
                   getN(R_dada),
                   getN(merged),
                   rowSums(nochim))
  }
  colnames(trackReads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
  trackReads <- trackReads/trackReads[,1]  # calculate relative values
  trackReads <- melt(trackReads)
  colnames(trackReads) <- c("name", "step", "reads")

  # save output
  saveRDS(nochim, file = seqtab_out)    # ASV table
  saveRDS(tax, file = taxtab_out)       # taxonomy table
  saveRDS(trackReads, file = track_out) # track reads
}
