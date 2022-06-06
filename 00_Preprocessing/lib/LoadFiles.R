LoadFiles <- function(data, file.type, fastq.processed){
    if(file.type == "BIGWIG"){
        df <- plyranges::read_bigwig(data)

    } else if(file.type == "BED" | file.type == "CSV" | file.type == "TXT"){
        if(fastq.processed){
            df <- read_bed(data)
        } else {
            df <- fread(data)

            if(any(c("Watson", "Crick") %in% colnames(df))){
                setnames(df, paste0("V", 1:4))
                df[, `:=`(V3 = (V2+1)-V2, V4 = ifelse(V3 == 1, "+", "-"))]
            } else if(length(grep(pattern = "Recombination", x = data)) > 0){
                setnames(df, paste0("V", 1:4))
                df <- df[, .(V1, V2)]
                df[, `:=`(V3 = (V2+1)-V2, V4 = "+")]
            } else {
                setnames(df, paste0("V", 1:ncol(df)))
                
                if(any(!is.na(str_extract(string = unique(df$V1), pattern = "chr")))){
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
            df <- as_granges(df)
        }
    } else if(file.type == "BEDGRAPH"){
        df <- plyranges::read_bed_graph(data)

    } else if(file.type == "NARROWPEAK"){
        df <- plyranges::read_narrowpeaks(data)
        
    }

    df <- df %>% 
        dplyr::filter(seqnames %in% paste0("chr", 1:22)) %>% 
        dplyr::group_by(seqnames)

    if(file.type != "BED" & file.type != "CSV" & file.type != "TXT"){
        df <- df %>% 
            dplyr::arrange(dplyr::desc(score))
    }
    
    return(df)
}