PreprocessFiles <- R6::R6Class(
    classname = "PreprocessFiles",
    public = list(
        initialize = function(){
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
                "hg18", "hg19", "hg38", "hs37d5", 
                "1000_Genomes_exp", "1000_Genomes_Pilot"
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
                    system(paste0("/usr/local/bin/wget https://ftp.ncbi.nlm.nih.gov/1000genomes/", 
                                  "ftp/technical/reference/human_g1k_v37.fasta.gz"))
                    system("/bin/mv human_g1k_v37.fasta.gz ../../data/ref/1000_Genomes_exp/")
                    genome <- Biostrings::readDNAStringSet(
                        filepath = "../../data/ref/1000_Genomes_exp/human_g1k_v37.fasta.gz"
                    )
                } else if("1000_Genomes_Pilot" == dir.to.create[i]){
                    cat(cur.msg, l, "\n", sep = "")
                    system(paste0("/usr/local/bin/wget https://ftp.ncbi.nlm.nih.gov/1000genomes/", 
                                  "ftp/pilot_data/technical/reference/human_b36_male.fa.gz"))
                    system(paste0("/usr/local/bin/wget https://ftp.ncbi.nlm.nih.gov/1000genomes/", 
                                  "ftp/pilot_data/technical/reference/human_b36_male.fa.gz.fai"))
                    system("/bin/mv human_b36_male.fa.gz* ../../data/ref/1000_Genomes_Pilot/")
                    genome.index <- Biostrings::fasta.index(
                        filepath = "../../data/ref/1000_Genomes_Pilot/human_b36_male.fa.gz"
                    )
                    genome <- Biostrings::readDNAStringSet(
                        filepath = "../../data/ref/1000_Genomes_Pilot/human_b36_male.fa.gz"
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
                } else if("1000_Genomes_Pilot" == dir.to.create[i]){
                    system("/bin/rm ../../data/ref/", 
                           "1000_Genomes_Pilot/human_b36_male.fa.gz")
                }

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }

            final.t <- Sys.time() - start.time
            cat("Final time taken:", final.t[[1]], attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Function that loads and preprocesses the files.
        #' @return None.
        load_and_process_files = function(){
            start.time <- Sys.time()

            for(i in 1:nrow(private$org_file)){
                # progress message
                t1 <- Sys.time()
                cur.msg <- paste0("Processing file ", i, " of ", nrow(private$org_file))
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))
                
                # obtain file name
                path.to.save.bp <- paste0(
                    "../../data/", private$org_file[i, `Fragmentation type`],
                    "/", private$org_file[i, `Experiment folder`]
                )
                all.files <- list.files(path = path.to.save.bp, full.names = TRUE)
                files.to.process <- all.files[!file.info(all.files)$isdir]
                files.to.process <- stringr::str_sort(files.to.process, numeric = TRUE)
                do.not.process <- stringr::str_extract(
                    string = files.to.process, 
                    pattern = "fastq|fasta|bam|bai"
                )
                files.to.process <- files.to.process[is.na(do.not.process)]
                file.names <- stringr::str_sort(stringr::str_extract(
                    string = files.to.process, 
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
                        file.type = private$org_file[i, `File format`],
                        fastq.processed = as.logical(private$org_file[i, `Fastq Processed`])
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

                    dir.create(paste0(path.to.save.bp, "/breakpoint_positions/"), showWarnings = FALSE)
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
                            strands = strands,
                            fastq.processed = as.logical(private$org_file[i, `Fastq Processed`])
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
                                strands = strands,
                                fastq.processed = as.logical(private$org_file[i, `Fastq Processed`])
                            )
                        }
                    }
                }

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }

            final.t <- Sys.time() - start.time
            cat("Final time taken:", final.t[[1]], attr(final.t, "units"), "\n")
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
            df <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )
            private$org_file <- df[Processed == "TRUE" & `DSB Map` == "TRUE"]
        },

        #' @description
        #' Load files for preprocessing.
        #' @param file Character vector of path to file.
        #' @param file.type Character vector of the file type.
        #' @param fastq.processed Boolean of whether the fastq had to be processed.
        #' @return Tibble.
        get_file = function(file, file.type, fastq.processed){
            if(file.type == "BIGWIG"){
                df <- plyranges::read_bigwig(file)

            } else if(file.type == "BED" | file.type == "CSV" | file.type == "TXT"){
                if(fastq.processed){
                    df <- plyranges::read_bed(file)
                } else {
                    df <- fread(file)

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
                }
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
        #' @param fastq.processed Boolean of whether the fastq had to be processed.
        #' @return None.
        process_file = function(data, chr, file.type, start.index, 
                                end.index, path.to.save.bp, both.ends, 
                                strands, fastq.processed){
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
                    if(fastq.processed){
                        dt <- data.table(
                            start.pos = ifelse(strand == "+",
                                start(filtered.ranges)+start.index,
                                end(filtered.ranges)+end.index)
                        )
                    } else {
                        dt <- data.table(
                            start.pos = ifelse(strand == "+",
                                start(filtered.ranges)+start.index,
                                start(filtered.ranges))
                        )
                    }
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
        }

    )
)