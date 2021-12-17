#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

g4_outfile <- args[1]
region_file <- args[2] ## col1 chrom, col2 start, col3 end
separate_by <- args[3]
first_position_ignore <- args[4]
out_file <- args[5]
output_include_origin_peak <- args[6]

# g4_outfile <- "/fh/fast/gottardo_r/yezheng_working/G-quadruplex/results/CSV/dm6_g4_grouped_window10_thres1.5.csv"
# region_file <- "/fh/fast/gottardo_r/yezheng_working/G-quadruplex/data/AllgenesChrStartEndNameStrand.csv"
# separate_by <- ","
# first_position_ignore <- FALSE
# out_file <- "/fh/fast/gottardo_r/yezheng_working/G-quadruplex/results/CSV/dm6_g4_grouped_window10_thres1.5_motif.bed"
# output_include_origin_peak <- FALSE

g4_out <- read.table(g4_outfile, header = T, sep = "\t")
region <- read.table(region_file, sep = separate_by, header = F, fill = T)

if(first_position_ignore){
    region$width <- region$V3 - region$V2
}else{
    region$width <- region$V3 - region$V2 + 1
}

get_g4_hunter_motif <- function(g4_res, region) {
    region_seq_len <- cumsum(data.frame(region)$width)
    motif_coord <- c()
    ind <- 1
    for (i in 1:nrow(g4_res)) {
        while (g4_res$POSITION[i] > region_seq_len[ind]) {
            ind <- ind + 1
        }
        if(ind > 1){
            position <- g4_res$POSITION[i] - region_seq_len[ind - 1]
        }else{
            position <- g4_res$POSITION[i]
        }

        motif_coord <- rbind(motif_coord, data.frame(
            chrom = region[ind, ]$V1,
            region_start = region[ind, ]$V2,
            region_end = region[ind, ]$V3,
            motif_start = region[ind, ]$V2 + position,
            motif_end = region[ind, ]$V2 + position + g4_res$LENGTH[i] - 1,
            motif_length = g4_res$LENGTH[i],
            motif_score = g4_res$SCORE[i],
            motif_sequence = g4_res$SEQUENCE[i],
            motif_strand = gsub("1", "+", gsub("-1", "-", sign(g4_res$SCORE[i])))
        ))
    }
    return(motif_coord)
}
res <- get_g4_hunter_motif(g4_out, region)

if(output_include_origin_peak){
    write.table(res, file = out_file, quote = F, col.names = T, row.names = F, sep = "\t")
}else{

    write.table(res[, c("chrom", "motif_start", "motif_end", "motif_length", "motif_score", "motif_sequence", "motif_strand")], file = out_file, quote = F, col.names = T, row.names = F, sep = "\t")
}
