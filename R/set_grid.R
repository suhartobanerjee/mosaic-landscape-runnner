library(data.table)
library(GenomicRanges)


read_sv_data <- function(path = "../data/TALL03-DEA5.strict.filtered.txt") {
    
    raw <- fread(path)
    # removing inf llr_to_ref
    sv_dt <- raw[llr_to_ref != Inf]
    setnames(sv_dt,
             c("start", "end"),
             c("start_loc", "end_loc")
    )

    sv_gro <- makeGRangesFromDataFrame(sv_dt,
                                       seqnames.field = "chrom",
                                       start.field = "start_loc",
                                       end.field = "end_loc",
                                       keep.extra.columns = T
    )

    return(list(
        sv_gro = sv_gro,
        sv_dt = sv_dt
    ))
}


read_ideo_data <- function(path = "../data/ideogram_scaffold.tsv") {
    
    ideo_dt <- fread(path)

    ideo_gro <- makeGRangesFromDataFrame(ideo_dt,
                                         seqnames.field = "chrom",
                                         start.field = "start_loc",
                                         end.field = "end_loc",
                                         keep.extra.columns = T
    )

    return(list(
        ideo_gro = ideo_gro,
        ideo_dt = ideo_dt
    ))
}


set_grid <- function(sv_path = NULL) {

    if(is.null(sv_path)) {
        ret_list <- read_sv_data()
    } else {
        ret_list <- read_sv_data(sv_path)
    }
    sv_gro <- ret_list[[1]]
    sv_dt <- ret_list[[2]]

    ret_list <- read_ideo_data()
    ideo_gro <- ret_list[[1]]
    ideo_dt <- ret_list[[2]]

    sv_on_ideo_dt <- rbind(ideo_dt, sv_dt, fill = T)
    sv_on_ideo_gro <- makeGRangesFromDataFrame(sv_on_ideo_dt,
                                               seqnames.field = "chrom",
                                               start.field = "start_loc",
                                               end.field = "end_loc",
                                               keep.extra.columns = T
    )

    disjoin_gro <- disjoin(sv_on_ideo_gro)
    disjoin_gro$id <- c(1:length(disjoin_gro))

    disjoin_dt <- as.data.table(disjoin_gro)
    setnames(disjoin_dt,
             c("seqnames", "start", "end"),
             c("chrom", "start_loc", "end_loc")
    )

    hits <- findOverlaps(sv_gro, disjoin_gro)
    sv_dt[queryHits(hits), id := disjoin_gro[subjectHits(hits)]$id]

    ljoin_dt <- sv_dt[disjoin_dt, on = .(id)]
    ljoin_dt[is.na(chrom), `:=`(chrom = i.chrom,
                                start_loc = i.start_loc,
                                end_loc = i.end_loc)]
    ljoin_dt[is.na(sv_call_name), `:=`(sv_call_name = "ref",
                                       cell = "all",
                                       llr_to_ref = 0)]
    ljoin_dt[, width := end_loc - start_loc]
    ljoin_dt <- ljoin_dt[width != 0]
    grid_dt <- ljoin_dt[, .(chrom,
                            start_loc,
                            end_loc,
                            width,
                            cell,
                            sv_call_name,
                            llr_to_ref,
                            id)]

    return(grid_dt)
}

