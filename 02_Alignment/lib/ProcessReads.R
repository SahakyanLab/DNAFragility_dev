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
                len.of.loop <- 1:nrow(private$org_file)
            } else {
                len.of.loop <- self$which_exp_ind
            }

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()
                cur.msg <- paste0("Processing experiment ", i, 
                                  " of ", nrow(private$org_file))
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(cur.msg, l, "\n", sep = "")

                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )
                private$file_format <- private$org_file[i, `File format`]
                private$alignment_strand <- private$org_file[i, Alignment_strand]

                if((private$file_format == "BAM") & (length(list.files(
                        path = paste0("../../data/", private$bp_exp),
                        pattern = ".fasta.gz")) == 0)){
                    # extract bam files
                    private$extract_bam_files()

                    # update org_file.csv file
                    private$org_file[i, `File format`] <- "FASTA"
                    fwrite(
                        private$org_file,
                        file = "../../data/org_file.csv"
                    )
                    private$get_org_file()
                    private$file_format <- private$org_file[i, `File format`]
                }

                reads.already.processed <- list.files(
                    path = paste0("../../data/", private$bp_exp,
                                  "/breakpoint_positions"),
                    pattern = "alignment_file_[[:digit:]].txt",
                    recursive = TRUE
                )

                if(length(reads.already.processed) == 0){
                    # loop over each chromosome
                    for(chr in 1:22){
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
                            index <- max(floor(fasta.lines/self$interval), 1)
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
                private$calc_avg_levdist()

                # concatenate breakpoints into single file
                private$concat_breakpoints()   
                
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
            df <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )
            private$org_file <- df[Processed == "FALSE" & `DSB Map` == "TRUE"]
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
                    
                    # load bam with region limits
                    bam.content <- Rsamtools::scanBam(
                        file = bamFile, 
                        param = param)[[1L]] 
                } else {
                    # load bam with region limits
                    bam.content <- Rsamtools::scanBam(file = bamFile)[[1L]] 
                }
                
                # extract sequence data
                sequences <- bam.content$seq
                names(sequences) <- paste0(bam.content$pos, "_", bam.content$strand)
                
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
            }

            # locate BAM file
            files <- list.files(
                path = paste0("../../data/", private$bp_exp),
                pattern = "*.bam"
            )
            files <- stringr::str_sort(files, numeric = TRUE)

            if(length(files) > 2){
                files <- list.files(
                    path = paste0("../../data/", private$bp_exp),
                    pattern = "*.bam$"
                )
                files <- stringr::str_sort(files, numeric = TRUE)
                
                for(i in 1:length(files)){
                    # progress message
                    t1 <- Sys.time()
                    cur.msg <- paste0("Processing BAM file for chr", i)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(paste0(cur.msg, l))

                    bamFile <- Rsamtools::BamFile(
                        paste0("../../data/", private$bp_exp, "/", files[i])
                    )
                    bam <- load_chr(
                        seqname = i, 
                        bamFile = bamFile, 
                        with.index = TRUE
                    )

                    total.time <- Sys.time() - t1
                    cat("DONE! --", signif(total.time[[1]], 2), 
                        attr(total.time, "units"), "\n")
                }
            } else {
                # check if index file exists
                check.ind.exist <- stringr::str_extract(
                    string = files, 
                    pattern = "\\.([[:alnum:]]+)$"
                )
                with.index <- ifelse(".bai" %in% check.ind.exist, TRUE, FALSE)
                bamFile <- Rsamtools::BamFile(
                    paste0("../../data/", private$bp_exp, "/", files[1])
                )

                for(i in 1:22){
                    # progress message
                    t1 <- Sys.time()
                    cur.msg <- paste0("Processing BAM file for chr", i)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(paste0(cur.msg, l))

                    bam <- load_chr(
                        seqname = i, 
                        bamFile = bamFile, 
                        with.index = with.index
                    )

                    total.time <- Sys.time() - t1
                    cat("DONE! --", signif(total.time[[1]], 2), 
                        attr(total.time, "units"), "\n")
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
                reads <- stringr::str_sub(
                    string = private$ref,
                    start = df$V2,
                    end = df$V3
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
            ref.seq <- stringr::str_sub(
                string = private$ref, 
                start = read.start.pos, 
                end = read.end.pos
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
                Mean <- sum(lev.dist.df$lev.dist*lev.dist.df$count)/
                        sum(lev.dist.df$count)
                SD <- sqrt(sum((lev.dist.df$lev.dist)**2*lev.dist.df$count)/
                          (sum(lev.dist.df$count)-1))
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
            cur.msg <- paste0("Concat breakpoints into chr-sep files for chr", i)
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
                    df <- df[lev.dist < mean(lev.dist.df$SD1)]
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

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")            
        },

        #' @description
        #' Formats the breakpoint csv files into kmertone-ready files.
        #' @return None.
        format_file_for_kmertone = function(){
            t1 <- Sys.time()
            cur.msg <- paste0("Formatting kmertone-ready files for chr", i)
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
                                    "/kmertone/chr", i, ".txt")
                )
            }

            total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")            
        }
    )
)