#############################################################################################
AlignReads <- function(
  chromosome.nr,
  ind,
  fasta.lines,
  reads.path,
  interval,
  BAM = FALSE,
  alignment.strands = "both"
  ){

  # check input classes
  if(!is.numeric(chromosome.nr)) stop("chromosome.nr needs to be a numeric value.")
  if(!is.numeric(ind)) stop("ind needs to be a numeric value.")
  if(!is.numeric(fasta.lines)) stop("fasta.lines needs to be a numeric value.")
  if(!is.character(reads.path)) stop("reads.path needs to be a character vector.")
  if(!is.numeric(interval)) stop("interval needs to be a numeric value.")
  if(!is.logical(BAM)) stop("BAM needs to be TRUE/FALSE logic.")
  if(!is.character(alignment.strands)) stop("alignment.strands needs to be a character vector.")

  if(ind == 0){
    reads <- fasta.index(paste0("../data/reads/", reads.path, "/chr", chromosome.nr, ".fasta.gz"),
                        skip = 0, 
                        nrec = interval)
  } else {
    # load read sequences
    to.skip <- seq(from = interval+1, to = fasta.lines, by = interval+1) 
    reads <- fasta.index(paste0("../data/reads/", reads.path, "/chr", chromosome.nr, ".fasta.gz"),
                        skip = to.skip[ind], 
                        nrec = interval)
  } 
  reads <- readDNAStringSet(reads)

  # extract length of read
  read.length <- width(reads)

  # extract position on chromosome
  read.names <- names(reads)

  if(BAM){
    # raw sequencing files in BAM format
    str.extract <- str_split(string = read.names, pattern = "_")

    # strand information
    strands <- sapply(str.extract, "[[", 2)

    if(alignment.strands == "plus"){
      reads <- reads[which(strands == "+")]
    } else if(alignment.strands == "minus"){
      reads <- reverseComplement(reads[which(strands == "-")])
    } else if(alignment.strands == "both"){
      to.revcomp <- which(strands == "-")
      reads <- c(reads[-to.revcomp], reverseComplement(reads[to.revcomp]))
    }

    # bp start position
    read.names <- names(reads)
    read.length <- width(reads)
    str.extract <- str_split(string = read.names, pattern = "_")
    read.start.pos <- sapply(str.extract, "[[", 1)
  } else {
    # raw sequencing files in fasta format
    read.start.pos <- str_extract(string = read.names, pattern = "(?<=pos=).*(?= mapq)")
  }

  read.start.pos <- as.integer(read.start.pos)
  read.end.pos <- read.start.pos+read.length-1

  if(!BAM | (alignment.strands == "both")){
    # get reverse complement of reads
    reads.revcomp <- reverseComplement(reads)
    reads.revcomp <- paste(reads.revcomp)
  }

  reads <- paste(reads)

  # get substrings of reference sequence
  ref.seq = str_sub(string = ref.seq.original, start = read.start.pos, end = read.end.pos)

  if(!BAM | (alignment.strands == "both")){
    outMatrix <- matrix(
      data = c(reads, reads.revcomp, ref.seq), 
      ncol = 3
    )
  } else {
    outMatrix <- matrix(
      data = c(reads, ref.seq), 
      ncol = 2
    )
  }

  # Levenshtein distance calculations
  lev.dist <- LevenshteinLoop(outMatrix)

  if(!BAM | (alignment.strands == "both")){
    best.align <- max.col(-lev.dist)

    # align first 2 nucleotides of read against reference sequence
    results <- ifelse(
      rep.int(best.align == 1, 2),
      ifelse(
        rep.int(
          str_sub(string = reads, start = 1, end = 2) == 
          str_sub(string = ref.seq, start = 1, end = 2), 
          times = 2
        ),
        c(read.start.pos-1, lev.dist[, 1]),
        NA
      ),
      ifelse(
        rep.int(
          str_sub(string = reads.revcomp, start = read.length-1, end = read.length) == 
          str_sub(string = ref.seq, start = read.length-1, end = read.length), 
          times = 2
        ),
        c(read.end.pos, lev.dist[, 2]),
        NA
      )
    )

    results <- na.omit(results)
    results.length <- length(results)
    split.results <- results.length/2

    df <- data.table(
      bp.start.pos = results[1:split.results],
      lev.dist = results[(split.results+1):results.length]
    )
  } else {
    df <- data.table(
      bp.start.pos = ifelse(alignment.strands == "plus", read.start.pos, read.end.pos),
      lev.dist = lev.dist
    )
    
    setnames(df, c("bp.start.pos", "lev.dist"))
  }

  df[, freq := .N, keyby = .(bp.start.pos, lev.dist)]

  df <- df[df[, .I[which.min(lev.dist)], by = bp.start.pos][["V1"]]]
  
  return(df)
}
#############################################################################################