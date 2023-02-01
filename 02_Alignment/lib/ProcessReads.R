ProcessReads <- R6::R6Class(
    classname = "ProcessReads",
    public = list(
        #' @field interval Numeric vector of how many reads to load into RAM.
        interval = 1000000,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        initialize = function(interval, which_exp_ind){
            if(!missing(interval)) self$interval <- interval
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Processes reads from each experiment
        #' @return None.
        process_reads = function(){
            if(is.null(self$which_exp_ind)){
                self$which_exp_ind <- which(
                    (private$org_file$Processed == "FALSE") & 
                    (private$org_file$`DSB Map` == "TRUE")
                )
            }
            len.of.loop <- self$which_exp_ind

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()

                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )		        
                private$file_format <- private$org_file[i, `File format`]
                private$alignment_strand <- private$org_file[i, Alignment_strand]

                cur.msg <- paste0("Processing exp ", i, "/", 
                                  nrow(private$org_file),
                                  ". ", private$bp_exp)
                cat(cur.msg, "\n", sep = "")

                if((private$file_format == "BAM") & (length(list.files(
                        path = paste0("../../data/", private$bp_exp),
                        pattern = ".fasta.gz")) == 0)){
                    # extract bam files
                    private$extract_bam_files()

                    # update org_file.csv file
                    private$org_file$`File format`[i] <- "FASTA"
                    fwrite(
                        private$org_file,
                        file = "../../data/org_file.csv"
                    )
                    private$get_org_file()
                    private$file_format <- private$org_file[i, `File format`]
                }
                
                # check if any alignments have been performed already
                raw.chr.alignments.done <- length(list.dirs(
                    path = paste0("../../data/", private$bp_exp,
                                "/breakpoint_positions")
                ))-1            
                concat.chr.alignments.done <- length(list.files(
                    path = paste0("../../data/", private$bp_exp,
                                  "/breakpoint_positions"),
                    pattern = "chr[0-9]|[10-22].csv",
                    recursive = TRUE
                ))              

                if(concat.chr.alignments.done < 22){
                    # loop over each chromosome
                    for(chr in max(c(raw.chr.alignments.done, 1)):22){
                        # progress message
                        t1 <- Sys.time()
                        cur.msg <- paste0("Aligning reads to reference for chr", chr)
                        l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                        cat(paste0(cur.msg, l))

                        private$exp_path <- paste0("../../data/", private$bp_exp, 
                                                   "/chr", chr, ".fasta.gz")
                        fasta.lines <- system(
                            command = paste0("/usr/bin/zgrep -Ec '>' ", private$exp_path),
                            intern = TRUE
                        )
                        fasta.lines <- as.numeric(fasta.lines)
                        if(is.null(self$interval)){
                            self$interval <- fasta.lines
                            index <- 0
                        } else {
                            index <- floor(fasta.lines/self$interval)
                        }
                        private$get_ref(ind = i, chr = chr)

                        for(id in 0:index){
                            private$align_reads(
                                chr = chr,
                                id = id,
                                fasta.lines = fasta.lines
                            )
                        }

                        total.time <- Sys.time() - t1
                        cat("DONE! --", signif(total.time[[1]], 2), 
                            attr(total.time, "units"), "\n")
                    }
                }

                # obtain average levdist for all chromosomes
                private$bp_pos_path <- paste0("../../data/", private$bp_exp,
                                              "/breakpoint_positions")
                if(raw.chr.alignments.done == 22) private$calc_avg_levdist()

                # concatenate breakpoints into single file
                if(concat.chr.alignments.done < 22) private$concat_breakpoints()

                # if ancient DNA, need to filter by the overlaps with the provided
                # bed files as per the publications. Else, exclude ENCODE blacklist regions.
                private$filter_dna_with_masks()
                
                # format files for kmertone
                private$format_file_for_kmertone()              

                # time taken for full processing for this experiment
                final.t <- Sys.time() - start.time
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
                cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                    attr(final.t, "units"), "\n")
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            }
        }
    ),
    private = list(
        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

        #' @field ref DNAStringSet of a human reference genome chromosome.
        ref = NULL,

        #' @field bp_exp Character vector of fragment type and experiment.
        bp_exp = NULL,

        #' @field bp_pos_path Character vector of path to breakpoint file.
        bp_pos_path = NULL,

        #' @field exp_path Character vector of path to raw seq fasta file.
        exp_path = NULL,

        #' @field file_format Character vector of the raw seq file format.
        file_format = NULL,

        #' @field alignment_strand Character vector of c("both", "plus", "minus").
        alignment_strand = NULL,

        #' @description
        #' Import full org_file.csv and filter for rows to be processed.
        #' @return None.
        get_org_file = function(){
            private$org_file <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )            
        },

        #' @description
        #' Get human reference genome file.
        #' @param ind Numeric vector of index subsetting org_file.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        get_ref = function(ind, chr){
            ref.seq <- private$org_file[ind, `Reference genome folder`]
            private$ref <- Biostrings::readDNAStringSet(
                filepath = paste0("../../data/ref/", ref.seq, 
                                  "/chr", chr, ".fasta.gz")
            )
        },

        #' @description
        #' Extract each chromosome from BAM files and save as 
        #' chromosome-separated fasta files.
        #' @return None.
        extract_bam_files = function(){
            load_chr <- function(seqname, bamFile, with.index){
                # reference names and lengths
                bamInfo <- seqinfo(bamFile)
                
                # start and end region to extract
                afrom <- 1 
                ato <- seqlengths(bamInfo[as.character(seqname)]) 
                
                if(with.index){
                    # set parameters for loading
                    param <- Rsamtools::ScanBamParam(
                        what = c('seq', 'pos', 'strand'),
                        which = GRanges(seqname, IRanges(afrom, ato)),
                        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
                    )
                    
                    # progress message
                    t1 <- Sys.time()
                    cur.msg <- "Loading in bam file"
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(paste0(cur.msg, l))     

                    # load bam with region limits
                    bam.content <- Rsamtools::scanBam(
                        file = bamFile, 
                        param = param)[[1L]]

                    total.time <- Sys.time() - t1
                    cat("DONE! --", signif(total.time[[1]], 2), 
                        attr(total.time, "units"), "\n")                    
                } else {
                    # load bam with region limits
                    bam.content <- Rsamtools::scanBam(file = bamFile)[[1L]] 
                }
                
                # extract sequence data
                sequences <- bam.content$seq
                names(sequences) <- paste0(bam.content$pos, "_", bam.content$strand)
                
                # progress message
                t1 <- Sys.time()
                cur.msg <- "Saving as fasta file"
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l)) 

                # save as fasta file
                dir.create(
                    path = paste0("../../data/", private$bp_exp),
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                Biostrings::writeXStringSet(
                    x = sequences, 
                    filepath = paste0("../../data/", private$bp_exp, 
                                      "/chr", as.character(seqname), 
                                      ".fasta.gz"), 
                    format = "fasta", 
                    compress = TRUE
                )

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")   
            }

            # locate BAM file
            files <- list.files(
                path = paste0("../../data/", private$bp_exp, "/bam"),
                pattern = "*.bam",
                full.names = TRUE
            )
            files <- stringr::str_sort(files, numeric = TRUE)

            if(length(files) > 2){
                files <- list.files(
                    path = paste0("../../data/", private$bp_exp, "/bam"),
                    pattern = "*.bam$",
                    full.names = TRUE
                )
                files <- stringr::str_sort(files, numeric = TRUE)
                
                for(i in 1:length(files)){
                    # progress message
                    t1 <- Sys.time()
                    cur.msg <- paste0("Processing BAM file for chr", i)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, "\n", sep = "")

                    bamFile <- Rsamtools::BamFile(files[i])
                    bam <- load_chr(
                        seqname = i, 
                        bamFile = bamFile, 
                        with.index = TRUE
                    )
                }
            } else {
                # check if index file exists
                check.ind.exist <- stringr::str_extract(
                    string = files, 
                    pattern = "\\.([[:alnum:]]+)$"
                )
                with.index <- ifelse(".bai" %in% check.ind.exist, TRUE, FALSE)
                bamFile <- Rsamtools::BamFile(files[1])

                for(i in 1:22){
                    # progress message
                    t1 <- Sys.time()
                    cur.msg <- paste0("Processing BAM file for chr", i)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, "\n", sep = "")   

                    bam <- load_chr(
                        seqname = i, 
                        bamFile = bamFile, 
                        with.index = with.index
                    )
                }
                # Cleanup
                # close(bamFile)
            }
        },

        #' @description
        #' Aligns sequencing reads to human reference genome. Finds best match.
        #' @param chr Numeric vector of chromosome number.
        #' @param id Numeric vector of number of subsets of raw seq read.
        #' @param fasta.lines Numeric vector of number of raw seq reads in chromosome.
        align_reads = function(chr, id, fasta.lines){
            if(private$file_format == "BED"){
                df <- fread(
                    paste0(private$exp_path, "/chr", chr, ".bed"),
                    showProgress = FALSE
                )
                df[, `:=`(V4 = NULL, V5 = NULL)]

                # get substrings of reference sequence
                reads <- substring(
                    text = private$ref,
                    first = df$V2,
                    last = df$V3
                )

                if(private$alignment_strand == "both"){
                    reads.revcomp <- Biostrings::reverseComplement(
                        Biostrings::DNAStringSet(reads)
                    )
                    reads.revcomp <- paste(reads.revcomp)
                    outMatrix <- matrix(
                        data = c(reads, reads.revcomp, ref.seq),
                        ncol = 3
                    )
                }
            }

            if(id == 0){
                reads <- Biostrings::fasta.index(
                    private$exp_path,
                    skip = 0, 
                    nrec = self$interval
                )
            } else {
                # load read sequences
                to.skip <- seq(
                    from = self$interval+1, 
                    to = fasta.lines, 
                    by = self$interval+1
                )
                reads <- Biostrings::fasta.index(
                    private$exp_path,
                    skip = to.skip[id], 
                    nrec = self$interval
                )
            } 
            reads <- Biostrings::readDNAStringSet(reads)

            # extract length of read
            read.length <- width(reads)

            # extract position on chromosome
            read.names <- names(reads)

            if(private$file_format == "BAM"){
                # raw sequencing files in BAM format
                str.extract <- stringr::str_split(
                    string = read.names, 
                    pattern = "_"
                )

                # strand information
                strands <- sapply(str.extract, "[[", 2)

                if(private$alignment_strand == "plus"){
                    reads <- reads[which(strands == "+")]
                } else if(private$alignment_strand == "minus"){
                    reads <- Biostrings::reverseComplement(
                        reads[which(strands == "-")]
                    )
                } else if(private$alignment_strand == "both"){
                    to.revcomp <- which(strands == "-")
                    reads <- c(
                        reads[-to.revcomp], 
                        Biostrings::reverseComplement(reads[to.revcomp])
                    )
                }

                # bp start position
                read.names <- names(reads)
                read.length <- width(reads)
                str.extract <- stringr::str_split(
                    string = read.names, 
                    pattern = "_"
                )
                read.start.pos <- sapply(str.extract, "[[", 1)
            } else {
                # raw sequencing files in fasta format
                # temp for when read.names equals
                # "60001"
                if(str_detect(read.names[1], "^[:digit:]+$")){
                    read.start.pos <- read.names
                } else {
                    # temp for when read.names equals
                    # "ERR1025627.55/2 primary ref=chr1 pos=10001 mapq=254"
                    temp <- stringr::str_detect(
                        string = read.names[1], 
                        pattern = "pos"
                    )
                    if(temp){
                        read.start.pos <- stringr::str_extract(
                            string = read.names, 
                            pattern = "(?<=pos=).*(?= mapq)"
                        )
                    } else {
                        # temp for when read.names equals
                        # "60005_+"
                        temp <- stringr::str_replace_all(
                            string = read.names[1], 
                            pattern = "[^[:digit:]]", 
                            replacement = ""
                        )
                        
                        if(str_detect(temp, "^[:digit:]+$")){
                            read.start.pos <- stringr::str_replace_all(
                                string = read.names, 
                                pattern = "[^[:digit:]]", 
                                replacement = ""
                            )
                        }
                    }
                }
            }
            read.start.pos <- as.integer(read.start.pos)
            read.end.pos <- read.start.pos+read.length-1

            if(private$file_format != "BAM" | (private$alignment_strand == "both")){
                # get reverse complement of reads
                reads.revcomp <- Biostrings::reverseComplement(reads)
                reads.revcomp <- paste(reads.revcomp)
            }
            reads <- paste(reads)

            # get substrings of reference sequence
            ref.seq <- substring(
                text = private$ref, 
                first = read.start.pos, 
                last = read.end.pos
            )

            if(private$file_format != "BAM" | (private$alignment_strand == "both")){
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

            if(private$file_format != "BAM" | (private$alignment_strand == "both")){
                best.align <- max.col(-lev.dist)

                # align first 2 nucleotides of read against reference sequence
                results <- ifelse(
                    rep.int(best.align == 1, 2),
                    ifelse(
                        rep.int(
                            stringr::str_sub(string = reads, 
                                             start = 1, 
                                             end = 2) == 
                            stringr::str_sub(string = ref.seq, 
                                             start = 1, 
                                             end = 2), 
                            times = 2
                        ),
                        c(read.start.pos-1, lev.dist[, 1]),
                        NA
                    ),
                    ifelse(
                        rep.int(
                            stringr::str_sub(string = reads.revcomp, 
                                             start = read.length-1, 
                                             end = read.length) == 
                            stringr::str_sub(string = ref.seq, 
                                             start = read.length-1, 
                                             end = read.length),
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
                    bp.start.pos = ifelse(private$alignment_strand == "plus", 
                                          read.start.pos, 
                                          read.end.pos),
                    lev.dist = lev.dist
                )
                setnames(df, c("bp.start.pos", "lev.dist"))
            }
            df[, freq := .N, keyby = .(bp.start.pos, lev.dist)]
            df <- df[df[, .I[which.min(lev.dist)], by = bp.start.pos][["V1"]]]
            
            dir.create(
                paste0("../../data/", private$bp_exp,
                       "/breakpoint_positions/chr", 
                       chr),
                recursive = TRUE,
                showWarnings = FALSE
            )
            fwrite(
                df,
                file = paste0(
                    "../../data/", private$bp_exp,
                    "/breakpoint_positions/chr", chr, 
                    ifelse(private$alignment_strand == "plus", "/plus_alignment_", 
                    ifelse(private$alignment_strand == "minus", "/minus_alignment_", 
                    "/alignment_")), 
                    "file_", id, ".txt"
                )
            )
        },

        #' @description
        #' Calcualtes the average levenstein distance per chromosome and plots results.
        #' @return None.
        calc_avg_levdist = function(){
            # progress bar
            t1 <- Sys.time()
            cur.msg <- "Calculating average levenshtein distance"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            lev.dist <- lapply(1:22, function(i){
                files <- list.files(
                    path = paste0(private$bp_pos_path, "/chr", i),
                    pattern = "alignment_file",
                    full.names = TRUE
                )
                files <- stringr::str_sort(files, numeric = TRUE)
                
                tables <- lapply(files, function(x){
                    fread(
                        file = x, 
                        header = TRUE,
                        drop = "bp.start.pos",
                        showProgress = FALSE
                    )
                })
                df <- rbindlist(tables)
                setnames(df, c("lev.dist", "freq"))
                          
                if(anyNA(df)){
                    df <- df[complete.cases(df)]
                }
                if("character" %in% sapply(df, class)){
                    df[, names(df) := lapply(.SD, as.numeric)]
                }

                # return frequency average and st.dev levenshtein distance
                lev.dist.df <- df[, .(count = .N), by = lev.dist]
                Mean <- sum(lev.dist.df$lev.dist*lev.dist.df$count, na.rm = TRUE)/
                        sum(lev.dist.df$count, na.rm = TRUE)
                SD <- sqrt(sum((lev.dist.df$lev.dist)**2*lev.dist.df$count, na.rm = TRUE)/
                          (sum(lev.dist.df$count, na.rm = TRUE)-1))
                return(list(Mean, SD))
            })

            # concatenate results and calculate overall statistics
            df <- rbindlist(lev.dist)

            if("list" %in% sapply(df, class)){
                df[, names(df) := lapply(.SD, as.numeric)]
            }
            setnames(df, c("Mean", "SD"))
            df[, `:=`(SD1 = Mean+SD, SD2 = Mean+2*SD)]

            # rename column
            df[, Chromosomes := 1:nrow(df)]
            setorder(df, -Chromosomes)
            df[, Chromosomes := paste0("Chr", 22:1)]
            df[, `:=`(Chromosomes = forcats::fct_inorder(Chromosomes))]

            # obtain bar plot of levenshtein distance
            lev.plot <- ggplot(df) +
                geom_bar(aes(
                    x = Chromosomes,
                    y = Mean),
                    stat = "identity",
                    fill = "skyblue",
                    alpha = 0.7) + 
                geom_pointrange(aes(
                    x = Chromosomes,
                    y = Mean,
                    ymin = Mean,
                    ymax = SD1),
                    color = "orange",
                    alpha = 1,
                    size = 0.8) +
                coord_flip() +
                labs(
                    title = "Average Levenshtein Distance plotted as Mean + 1 St.Dev", 
                    subtitle = paste0("Overall Average: ", signif(mean(df$Mean), 3))
                )

            dir.create(
                paste0("../figures/", private$bp_exp),
                showWarnings = FALSE,
                recursive = TRUE
            )
            ggsave(
                filename = paste0("../figures/", private$bp_exp, 
                                  "/AvgLevenshteinDistance.pdf"),
                plot = lev.plot,
                height = 8,
                width = 7
            )

            # save lev dist 
            dir.create(
                paste0("../../data/", private$bp_exp, "/average_levdist"),
                showWarnings = FALSE,
                recursive = TRUE
            )
            fwrite(
                x = df, 
                file = paste0("../../data/", private$bp_exp, 
                              "/average_levdist/AvgLevenshteinDistance.csv"), 
                row.names = FALSE
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Avg levenshtein distance for chr1 - chr22:", 
                signif(mean(df$Mean, na.rm = TRUE), digits = 3), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Concatenates all breakpoint files "alignment_file_i.txt" into 
        #' chromosome-separated csv files.
        #' @return None.
        concat_breakpoints = function(){
            t1 <- Sys.time()
            cur.msg <- "Concat breakpoints into chr-sep files"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            for(i in 1:22){
                # load data sets
                files <- list.files(
                    path = paste0(private$bp_pos_path, "/chr", i),
                    pattern = "alignment_file",
                    full.names = TRUE
                )
                files <- stringr::str_sort(files, numeric = TRUE)
                
                if(length(files) > 0){
                    # concatenate files
                    df <- lapply(files, function(x){
                        fread(
                            file = x, 
                            header = TRUE, 
                            showProgress = FALSE
                        )
                    })
                    df <- rbindlist(df)
                    setnames(df, c("start.pos", "lev.dist", "freq"))
                    df <- unique(df, by = "start.pos")
                    
                    # save bottom ~95% of the levenshtein distance score
                    lev.dist.file <- paste0("../../data/", private$bp_exp, 
                                            "/average_levdist/AvgLevenshteinDistance.csv")
                    if(file.exists(lev.dist.file)){
                        lev.dist.df <- fread(file = lev.dist.file, header = TRUE)
                        df <- df[lev.dist < mean(lev.dist.df$SD1, na.rm = TRUE)]
                    } else {
                        df <- df[lev.dist < 1]
                    }

                    fwrite(
                        df,
                        file = paste0("../../data/", private$bp_exp, 
                                        "/breakpoint_positions/chr", i, ".csv"),
                        row.names = FALSE
                    )
                    
                    # rm txt files and only keep csv files
                    folder.to.rm <- unique(stringr::str_remove(
                        string = files,
                        pattern = "/alignment_file_[[:digit:]]+.txt"
                    ))
                    system(paste0("/bin/rm -r ", folder.to.rm))
                }
            }   

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")            
        },

        #' @description
        #' If ancient DNA, need to filter by the overlaps with the provided
        #' bed files as per the publications. Else, exclude ENCODE blacklist regions.
        #' @return None.
        filter_dna_with_masks = function(){
            t1 <- Sys.time()
            cur.msg <- "Filter genomic DNA by provided masked/blacklist bed files"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # load data sets
            files <- list.files(
                path = private$bp_pos_path,
                pattern = "chr[0-9]|[10-22].csv",
                full.names = TRUE
            )
            files <- stringr::str_sort(files, numeric = TRUE)            

            if(length(files) > 0){
                for(i in 1:22){
                    df <- fread(files[i], showProgress = FALSE)

                    # regions of interest here
                    if(grepl(pattern = "Ancient_DNA", x = private$bp_exp)){
                        df.bp <- data.table(
                            seqnames = i, 
                            start = df$start.pos, 
                            width = 1
                        )
                        df.bp <- as_granges(df.bp)

                        # authors said these are regions to INCLUDE
                        df.bed <- plyranges::read_bed(paste0(
                            "../../data/", private$bp_exp, 
                            "/filterbed/chr", i, "_mask.bed"
                        ))
                        df.bed.all.pos <- plyranges::as_granges(df.bed)
                        # filter by overlaps
                        df.overlaps <- plyranges::filter_by_overlaps(df.bp, df.bed.all.pos)                        
                    } else {
                        # these are regions to EXCLUDE 
                        if(grepl(pattern = "hg38", 
                                 x = private$org_file$`Reference genome folder`[i], 
                                 ignore.case = TRUE)){
                            df.bed.all.pos <- plyranges::read_bed(
                                "../../data/blacklists/hg38.bed"
                            )                            
                        } else {
                            df.bed.all.pos <- plyranges::read_bed(
                                "../../data/blacklists/hg19.bed"
                            )
                        }
                        seqnames.convention <- any(grepl(
                            pattern = "chr", 
                            x = attr(attr(df.bed.all.pos, "seqnames"), "values")
                        ))
                        
                        # import true breakpoints
                        df.bp <- data.table(
                            seqnames = ifelse(seqnames.convention, paste0("chr", i), i), 
                            start = df$start.pos, 
                            width = 1
                        )
                        df.bp <- as_granges(df.bp)
                        df.overlaps <- plyranges::filter_by_non_overlaps(df.bp, df.bed.all.pos)
                    }

                    # as data.table
                    dt.overlaps <- as.data.table(df.overlaps)
                    dt.overlaps[, `:=`(
                        seqnames = NULL,
                        end = NULL, 
                        width = NULL, 
                        strand = NULL
                    )]
                    dt.overlaps <- as.data.table(unique(dt.overlaps$start))
                    setnames(dt.overlaps, "start.pos")

                    fwrite(
                        dt.overlaps,
                        file = files[i]
                    )
                }
            }
            
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")             
        },     

        #' @description
        #' Formats the breakpoint csv files into kmertone-ready files.
        #' @return None.
        format_file_for_kmertone = function(){
            t1 <- Sys.time()
            cur.msg <- "Formatting kmertone-ready files"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))
           
            for(i in 1:22){
                # obtain breakpoints
                fetch.file <- paste0(
                    "../../data/", private$bp_exp,
                    "/breakpoint_positions/chr", 
                    i, ".csv"
                )
                
                df <- fread(
                    file = fetch.file,
                    select = "start.pos",
                    showProgress = FALSE
                )
                df[, chromosome := rep(paste0("chr", i), nrow(df))]
                setcolorder(df, c("chromosome", "start.pos"))
                
                dir.create(
                    paste0("../../data/", private$bp_exp, "/kmertone"),
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                fwrite(
                    df, 
                    row.names = FALSE, 
                    file = paste0("../../data/", private$bp_exp, 
                                    "/kmertone/chr", i, ".csv")
                )
            }

            total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")            
        }
    )
)
