PreprocessFiles <- R6::R6Class(
    classname = "PreprocessFiles",
    public = list(
        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        initialize = function(which_exp_ind){
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Download packages or files of the human reference genome.
        #' @return None.
        get_genomes = function(){
            start.time <- Sys.time()
            packages <- c(
                "BSgenome.Hsapiens.UCSC.hg18",
                "BSgenome.Hsapiens.UCSC.hg19",
                "BSgenome.Hsapiens.UCSC.hg38",
                "BSgenome.Hsapiens.1000genomes.hs37d5"
            )

            dir.to.create <- c(
                "hg18", "hg19", "hg38", "hs37d5", "1000_Genomes_exp"
            )

            for(package in packages){
                if(!package %in% rownames(installed.packages())){
                    BiocManager::install(package)
                }
            }

            for(i in 1:length(dir.to.create)){
                # progress message
                t1 <- Sys.time()
                cur.msg <- paste0("Downloading human reference genome: ", dir.to.create[i])
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                dir.create(
                    path = paste0("../../data/ref/", dir.to.create[i], "/"),
                    recursive = TRUE,
                    showWarnings = FALSE
                )

                if("1000_Genomes_exp" == dir.to.create[i]){
                    cat(cur.msg, l, "\n", sep = "")
                    system(paste0("/usr/local/bin/wget ", 
                                  "https://ftp.ncbi.nlm.nih.gov/1000genomes/", 
                                  "ftp/technical/reference/human_g1k_v37.fasta.gz"))
                    system("/bin/mv human_g1k_v37.fasta.gz ../../data/ref/1000_Genomes_exp/")
                    genome <- Biostrings::readDNAStringSet(
                        filepath = "../../data/ref/1000_Genomes_exp/human_g1k_v37.fasta.gz"
                    )
                } else {
                    suppressPackageStartupMessages(suppressWarnings(
                        library(packages[i], character.only = TRUE)
                    ))
                }

                # extract all pairs of chromosomes
                if(i == 1){
                    genome <- BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18
                } else if(i == 2){
                    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
                } else if(i == 3){
                    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
                } else if(i == 4){
                    genome <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
                }

                check.files.exist <- list.files(
                    path = paste0("../../data/ref/", dir.to.create[i]),
                    pattern = "*.fasta.gz"
                )

                if(length(check.files.exist) <= 22){
                    for(x in 1:22){
                        if(("1000_Genomes_exp" == dir.to.create[i]) | 
                           ("1000_Genomes_Pilot" == dir.to.create[i])){
                            chr <- genome[i]
                        } else {
                            chr <- Biostrings::DNAStringSet(genome[[x]])
                        }

                        # save as fasta file
                        Biostrings::writeXStringSet(
                            chr,
                            filepath = paste0("../../data/ref/", dir.to.create[i], 
                                              "/chr", x, ".fasta.gz"), 
                            format = "fasta",
                            compress = TRUE
                        )
                    }
                }
                
                if("1000_Genomes_exp" == dir.to.create[i]){
                    system("/bin/rm ../../data/ref/", 
                           "1000_Genomes_exp/human_g1k_v37.fasta.gz")
                } 
               
                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }

            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 2), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },
        
        #' @description
        #' Function that loads and preprocesses the files.
        #' @return None.
        load_and_process_files = function(){
            start.time <- Sys.time()            
            if(is.null(self$which_exp_ind)){
                self$which_exp_ind <- which(
                    private$org_file$Processed == "TRUE"
                )
            }
            len.of.loop <- self$which_exp_ind

            for(i in len.of.loop){
                # progress message
                t1 <- Sys.time()
                bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )

                cur.msg <- paste0("Processing exp ", i, "/", 
                                  nrow(private$org_file),
                                  ". ", bp_exp)
                cat(cur.msg, "\n", sep = "")

                if(private$org_file$`RMSD?`[i]){
                    # obtain file name
                    path.to.save.bp <- paste0(
                        "../../data/", private$org_file[i, `Fragmentation type`],
                        "/", private$org_file[i, `Experiment folder`],
                        ifelse(private$org_file[i, `File format`] == "FASTA", "/filterbed", "")
                    )
                    all.files <- list.files(path = path.to.save.bp, full.names = TRUE)
                    files.to.process <- all.files[!file.info(all.files)$isdir]
                    files.to.process <- stringr::str_sort(files.to.process, numeric = TRUE)                  
                    do.not.process <- stringr::str_extract(
                        string = files.to.process, 
                        pattern = "fastq|fasta|bam|bai|sh"
                    )                    
                    files.to.process <- files.to.process[is.na(do.not.process)]
                    file.names <- stringr::str_sort(stringr::str_extract(
                        string = basename(files.to.process), 
                        pattern = "chr[[:alnum:]]+"), 
                        numeric = TRUE
                    )
                    both.ends <- ifelse(
                        private$org_file[i, Alignment_strand] == "both", 
                        TRUE, FALSE
                    )

                    for(x in 1:length(files.to.process)){
                        df <- private$get_file(
                            file = files.to.process[x],
                            file.type = private$org_file[i, `File format`]
                        )

                        if(all(grepl(pattern = "plus|minus", x = basename(files.to.process)))){
                            if(grepl(pattern = "plus", x = basename(files.to.process)[i])){
                                strands <- "plus"
                            } else if(grepl(pattern = "minus", x = basename(files.to.process)[i])){
                                strands <- "minus"
                            }
                        } else {
                            strands <- ""
                        }

                        dir.create(
                            path = paste0(path.to.save.bp, "/breakpoint_positions/"), 
                            showWarnings = FALSE
                        )
                        if(all(paste0("chr", 1:22) %in% file.names)){
                            private$process_file(
                                data = df, 
                                chr = stringr::str_extract(
                                    string = file.names[i], 
                                    pattern = "[[:digit:]]+"
                                ),
                                file.type = private$org_file[i, `File format`],
                                start.index = private$org_file[i, Start],
                                end.index = private$org_file[i, End],
                                path.to.save.bp = path.to.save.bp,
                                both.ends = both.ends,
                                strands = strands
                            )
                        } else {
                            for(chr in 1:22){
                                private$process_file(
                                    data = df,
                                    chr = chr,
                                    file.type = private$org_file[i, `File format`],
                                    start.index = private$org_file[i, Start],
                                    end.index = private$org_file[i, End],
                                    path.to.save.bp = path.to.save.bp,
                                    both.ends = both.ends,
                                    strands = strands
                                )
                            }
                        }
                    }

                    for(chr in 1:22){
                        dir.create(
                            path = paste0(
                                path.to.save.bp, 
                                "/kmertone/"
                            ), 
                            showWarnings = FALSE
                        )
                        if(strands != ""){
                            dt.plus <- fread(paste0(
                                path.to.save.bp, 
                                "/breakpoint_positions/plus_chr", 
                                chr, ".csv"
                            ), showProgress = FALSE)
                            dt.minus <- fread(paste0(
                                path.to.save.bp, 
                                "/breakpoint_positions/minus_chr", 
                                chr, ".csv"
                            ), showProgress = FALSE)
                            dt <- rbind(dt.plus, dt.minus)

                            if(!"strand" %in% colnames(dt)){
                                dt[, strand := "+"]    
                            }                        
                            setorder(dt, strand)
                            fwrite(
                                dt, 
                                file = paste0(
                                    path.to.save.bp, 
                                    "/breakpoint_positions/chr", 
                                    chr, ".csv"
                            ))
                            system(paste0(
                                path.to.save.bp, 
                                "/breakpoint_positions/plus_chr", 
                                chr, ".csv"
                            ))
                            system(paste0(
                                path.to.save.bp, 
                                "/breakpoint_positions/minus_chr", 
                                chr, ".csv"
                            ))
                        }

                        # save breakpoint positions for kmertone
                        dt <- fread(paste0(
                            path.to.save.bp, 
                            "/breakpoint_positions/chr", 
                            chr, ".csv"
                        ), showProgress = FALSE)
                        dt[, chromosome := paste0("chr", chr)]

                        if(!"strand" %in% colnames(dt)){
                            dt[, strand := "+"]    
                        }
                        setcolorder(dt, c("chromosome", "start.pos", "strand"))
                        fwrite(
                            dt, 
                            file = paste0(
                                path.to.save.bp, 
                                "/kmertone/chr", 
                                chr, ".csv"
                        ))
                    }
                } else {
                    private$separate_file_by_chr_no_rmsd(i = i)
                }
            }

            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        }
    ),
    private = list(
        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

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
        #' Load files for preprocessing.
        #' @param file Character vector of path to file.
        #' @param file.type Character vector of the file type.
        #' @return Tibble.
        get_file = function(file, file.type){
            if(file.type == "BIGWIG"){
                df <- plyranges::read_bigwig(file)

            } else if(file.type == "BED" | file.type == "CSV" | file.type == "TXT"){
                df <- fread(file, showProgress = FALSE)
                if(any(c("Watson", "Crick") %in% colnames(df))){
                    setnames(df, paste0("V", 1:4))
                    df[, `:=`(V3 = (V2+1)-V2, V4 = ifelse(V3 == 1, "+", "-"))]
                } else if(length(grep(pattern = "Recombination", x = file)) > 0){
                    setnames(df, paste0("V", 1:4))
                    df <- df[, .(V1, V2)]
                    df[, `:=`(V3 = (V2+1)-V2, V4 = "+")]
                } else {
                    setnames(df, paste0("V", 1:ncol(df)))
                    
                    if(any(!is.na(stringr::str_extract(string = unique(df$V1), 
                                                        pattern = "chr")))){
                        df <- df[V1 %in% paste0("chr", 1:22)]
                    } else {
                        df <- df[V1 %in% as.character(1:22)]
                        df[, V1 := paste0("chr", V1)]
                    }

                    df <- df[, .(V1, V2, V3)]
                    df[, `:=`(V4 = V3-V2, V5 = "+")]
                    df[, V3 := NULL]
                }

                setnames(df, c("seqnames", "start", "width", "strand"))
                df <- plyranges::as_granges(df)                
            } else if(file.type == "BEDGRAPH"){
                df <- plyranges::read_bed_graph(file)

            } else if(file.type == "NARROWPEAK"){
                df <- plyranges::read_narrowpeaks(file)
                
            }

            df <- df %>% 
                dplyr::filter(seqnames %in% paste0("chr", 1:22)) %>% 
                dplyr::group_by(seqnames)

            if(file.type != "BED" & file.type != "CSV" & file.type != "TXT"){
                df <- df %>% 
                    dplyr::arrange(dplyr::desc(score))
            }
            
            return(df)
        },

        #' @description
        #' Splits the file up into chromosome-separated files.
        #' @param data Tibble of file containing breakpoints for all chromosomes.
        #' @param chr Character vector of the chromosome
        #' @param file.type Character vector of the file type.
        #' @param start.index Numeric vector of the start index of a strand. 
        #' Some files need to shift the index by a position.
        #' @param end.index Numeric vector of the end index of a strand. 
        #' Some files need to shift the index by a position.
        #' @param path.to.save.bp Character vector of the path to the breakpoint file.
        #' @param both.ends Boolean. If TRUE, both strands processed, else just PLUS strand.
        #' @param strands Character vector of c("plus", "minus", "") strands.
        #' @return None.
        process_file = function(data, chr, file.type, start.index, 
                                end.index, path.to.save.bp, both.ends, strands){
            df.filtered <- data %>% 
                dplyr::filter(seqnames == paste0("chr", as.numeric(chr)))
            filtered.ranges <- ranges(df.filtered)

            if(strands == "minus" || strands == "plus"){
                dt.start <- data.table(
                    # start pos
                    start.pos = start(filtered.ranges)+start.index
                )

                dt.start <- unique(dt.start, by = "start.pos")
                dt.start[, `:=`(strand = ifelse(strands == "minus", "-", "+"))]

                if(both.ends){
                    dt.end <- data.table(
                        # end position
                        start.pos = end(filtered.ranges)+end.index
                    )

                    dt.end <- unique(dt.end, by = "start.pos")
                    dt.end[, `:=`(strand = ifelse(strands == "minus", "-", "+"))]

                    dt <- rbind(dt.start, dt.end)
                    setorder(dt, start.pos)
                } else {
                    dt <- dt.start
                    setorder(dt, start.pos)
                }
            } else {
                strand <- as.character(strand(df.filtered))
                if(all(c("+", "-") %in% unique(strand))){
                    dt <- data.table(
                        start.pos = ifelse(strand == "+",
                            start(filtered.ranges)+start.index,
                            start(filtered.ranges))
                    )
                } else {
                    dt.start <- data.table(
                        # start pos
                        start.pos = start(filtered.ranges)+start.index
                    )

                    if(both.ends){
                        dt.end <- data.table(
                            # end position
                            start.pos = end(filtered.ranges)+end.index
                        )
                        dt.end <- unique(dt.end, by = "start.pos")
                        dt <- rbind(dt.start, dt.end)
                        setorder(dt, start.pos)
                    } else {
                        dt <- dt.start
                        setorder(dt, start.pos)
                    }
                }
            }

            # in case any duplicates still present, this is a safety check
            dt <- as.data.table(unique(dt$start.pos))
            setnames(dt, "start.pos")

            if(strands == ""){
                fwrite(
                    dt, 
                    paste0(path.to.save.bp, "/breakpoint_positions/chr", 
                           chr, ".csv")
                )
            } else {
                fwrite(
                    dt, 
                    paste0(path.to.save.bp, "/breakpoint_positions/", 
                           strands, "_chr", chr, ".csv")
                )
            }
        },

        #' @description
        #' For files where no RMSD calculations will be performed, just
        #' split the files by chromosomes here.
        #' @return None.
        separate_file_by_chr_no_rmsd = function(i){
            files <- list.files(
                path = paste0("../../data/", private$org_file$`Fragmentation type`[i],
                              "/", private$org_file$`Experiment folder`[i]),
                pattern = "bed|narrowpeak|broadpeak|csv",
                full.names = TRUE,
                ignore.case = TRUE
            )
            pos <- regexpr("\\.([[:alnum:]]+)$", files)
            file.ext <- ifelse(pos > -1, 
                substring(text = files, first = pos+1), 
                ""
            )

            dir.create(
                path = paste0("../../data/", private$org_file$`Fragmentation type`[i],
                              "/", private$org_file$`Experiment folder`[i], 
                              "/breakpoint_positions"),
                showWarnings = FALSE
            )
            dir.create(
                path = paste0("../../data/", private$org_file$`Fragmentation type`[i],
                              "/", private$org_file$`Experiment folder`[i], 
                              "/kmertone"),
                showWarnings = FALSE
            )

            if(grepl(pattern = "csv", x = file.ext, ignore.case = TRUE)){
                df <- fread(files, showProgress = FALSE)
                setnames(df, c("V1", "V2", "V3"))
                df <- df[, .(V1, V2, V3)]
                setnames(df, c("Chromosome", "Start", "End"))

                # separate out the files by chromosome
                df.split <- split(df, by = "Chromosome")

                # only retain autosomes
                df.split <- df.split[grepl(
                    pattern = "chr[0-9]|[10-22]", 
                    x = names(df.split), 
                    ignore.case = TRUE
                )]

                res <- lapply(1:length(df.split), function(x){
                    fwrite(
                        df.split[[x]],
                        file = paste0(
                            "../../data/", private$org_file$`Fragmentation type`[i],
                            "/", private$org_file$`Experiment folder`[i], 
                            "/breakpoint_positions/", names(df.split)[[x]], ".csv")
                    )
                    fwrite(
                        df.split[[x]],
                        file = paste0(
                            "../../data/", private$org_file$`Fragmentation type`[i],
                            "/", private$org_file$`Experiment folder`[i], 
                            "/kmertone/", names(df.split)[[x]], ".csv")
                    )
                })
            } else if(grepl(pattern = "narrowpeak|bed|broadpeak", x = file.ext, ignore.case = TRUE)){
                if(grepl(pattern = "narrowpeak", x = file.ext, ignore.case = TRUE)){
                    df <- plyranges::read_narrowpeaks(files)
                } else if(grepl(pattern = "bed|broadpeak", x = file.ext, ignore.case = TRUE)){
                    df <- tryCatch({
                        plyranges::read_bed(files)
                    }, error = function(e){
                        temp <- fread(files, showProgress = FALSE)
                        temp <- temp[, .(V1, V2, V3)]
                        setnames(temp, c("seqnames", "start", "end"))
                        temp[, start := start+1] # bed-files are 0-indexed
                        return(plyranges::as_granges(temp))
                    })
                }
                if(ncol(attr(df, "elementMetadata")) > 0){
                    df <- df %>% dplyr::select(-colnames(attr(df, "elementMetadata")))
                }
                df <- df %>% dplyr::arrange(seqnames)

                # all.chr.names <- attr(attr(df, "seqnames"), "values")
                # all.chr.present <- grepl(
                #     pattern = "chr[0-9]|[10-22]", 
                #     x = all.chr.names, 
                #     ignore.case = TRUE
                # )

                res <- lapply(1:22, function(x){
                    temp <- df %>% 
                        dplyr::filter(seqnames == paste0("chr", x)) %>% 
                        dplyr::arrange(start)
                    dt <- as.data.table(temp)
                    colnames(dt) <- paste0("V", 1:ncol(dt))
                    dt <- dt[, .(V1, V2, V3)]
                    setnames(dt, c("Chromosome", "Start", "End"))
                    fwrite(
                        dt, 
                        file = paste0(
                            "../../data/", private$org_file$`Fragmentation type`[i],
                            "/", private$org_file$`Experiment folder`[i], 
                            "/breakpoint_positions/", paste0("chr", x), ".csv")
                    )                            
                    fwrite(
                        dt, 
                        file = paste0(
                            "../../data/", private$org_file$`Fragmentation type`[i],
                            "/", private$org_file$`Experiment folder`[i], 
                            "/kmertone/", paste0("chr", x), ".csv")
                    )
                })
            }  else {
                stop("Need to implement method.")
            }
        }
    )
)
