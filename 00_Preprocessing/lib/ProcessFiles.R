ProcessFiles <- function(data, chr, file.type, start.index, end.index, 
                         path.to.save.bp, both.ends, strands, fastq.processed){
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
        fwrite(dt, paste0(path.to.save.bp, "/breakpoint_positions/chr", chr, ".csv"))
    } else {
        fwrite(dt, paste0(path.to.save.bp, "/breakpoint_positions/", strands, "_chr", chr, ".csv"))
    }
}