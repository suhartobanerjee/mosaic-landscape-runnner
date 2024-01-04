library(data.table)
library(GenomicRanges)


set_dt_names <- function(gro_dt) {
    setnames(gro_dt,
        c("seqnames", "start", "end"),
        c("chrom", "start_loc", "end_loc"),
        skip_absent = T
    )


    return(gro_dt)
}


read_sv_data <- function(path = "../data/TALL03-DEA5.strict.filtered.txt") {
    raw <- fread(path)
    # removing inf llr_to_ref
    # removing complex svs
    sv_dt <- raw[llr_to_ref != Inf]
    sv_dt <- sv_dt[!str_detect(sv_call_name, regex("complex|idup"))]
    set_dt_names(sv_dt)
    sv_dt[, width := end_loc - start_loc]

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


# set_sv_dt <- function(sv_path = NULL) {
#
#     if(is.null(sv_path)) {
#         ret_list <- read_sv_data()
#     } else {
#         ret_list <- read_sv_data(sv_path)
#     }
#     sv_gro <- ret_list[[1]]
#     sv_dt <- ret_list[[2]]
#
#     ret_list <- read_ideo_data()
#     ideo_gro <- ret_list[[1]]
#     ideo_dt <- ret_list[[2]]
#
#     sv_on_ideo_dt <- rbind(ideo_dt, sv_dt, fill = T)
#     sv_on_ideo_gro <- makeGRangesFromDataFrame(sv_on_ideo_dt,
#                                                seqnames.field = "chrom",
#                                                start.field = "start_loc",
#                                                end.field = "end_loc",
#                                                keep.extra.columns = T
#     )
#
#     disjoin_gro <- disjoin(sv_on_ideo_gro)
#     disjoin_gro$id <- c(1:length(disjoin_gro))
#
#     disjoin_dt <- as.data.table(disjoin_gro)
#     setnames(disjoin_dt,
#              c("seqnames", "start", "end"),
#              c("chrom", "start_loc", "end_loc")
#     )
#
#     hits <- findOverlaps(sv_gro, disjoin_gro)
#     sv_dt[queryHits(hits), id := disjoin_gro[subjectHits(hits)]$id]
#
#     ljoin_dt <- sv_dt[disjoin_dt, on = .(id)]
#     ljoin_dt[is.na(chrom), `:=`(chrom = i.chrom,
#                                 start_loc = i.start_loc,
#                                 end_loc = i.end_loc)]
#     ljoin_dt[is.na(sv_call_name), `:=`(sv_call_name = "ref",
#                                        cell = "all",
#                                        llr_to_ref = 0)]
#     ljoin_dt[, width := end_loc - start_loc]
#     ljoin_dt <- ljoin_dt[width != 0]
#     sv_dt <- ljoin_dt[, .(chrom,
#                             start_loc,
#                             end_loc,
#                             width,
#                             cell,
#                             sv_call_name,
#                             llr_to_ref,
#                             id)]
#
#     landscape_obj <- new("Landscape",
#                          sv_dt = sv_dt,
#                          block = min(sv_dt$width),
#                          grid = data.table()
#     )
#
#     return(landscape_obj)
# }


bin_genome <- function(block) {
    ret_list <- read_ideo_data()
    ideo_gro <- ret_list[[1]]

    chrom_lengths <- width(ideo_gro)
    names(chrom_lengths) <- seqnames(ideo_gro)
    chrom_lengths

    # binning genome
    bin_gen_gro <- tileGenome(
        chrom_lengths,
        tilewidth = block,
        cut.last.tile.in.chrom = T
    )
    bin_gen_gro$bin_id <- c(1:length(bin_gen_gro))
    bin_gen_dt <- as.data.table(bin_gen_gro)
    # so that the end does not overlap with a sv
    # this prevents the edge case where the end is
    # a round number and it overlaps with a sv start
    # which includes the previous bin to overlap with the sv
    # even though it is a 1 base overlap.
    bin_gen_dt[, end := end - 1]

    # making the gro again with the above change
    bin_gen_gro <- makeGRangesFromDataFrame(bin_gen_dt)


    return(list(
        bin_gen_gro = bin_gen_gro,
        bin_gen_dt = bin_gen_dt
    ))
}


set_grid <- function(sv_path = NULL,
                     block = 1e5) {
    if (is.null(sv_path)) {
        ret_list <- read_sv_data()
    } else {
        ret_list <- read_sv_data(sv_path)
    }
    sv_gro <- ret_list[[1]]
    sv_dt <- ret_list[[2]]

    # binning the genome
    ret_list <- bin_genome(block)
    bin_gen_gro <- ret_list[[1]]
    bin_gen_dt <- ret_list[[2]]

    # finding overlaps between binned genome
    # and the sv_dt
    hits <- findOverlaps(sv_gro, bin_gen_gro)

    # adding bin_id col to sv_dt by cbinding
    cdt <- cbind(
        bin_gen_dt[subjectHits(hits), .(bin_id)],
        sv_dt[queryHits(hits)]
    )

    # right outer join to include all bin_gen_dt with
    # cdt(which is sv_dt + bin_id)
    grid_dt <- cdt[bin_gen_dt, on = .(bin_id)]

    # unique names for the bin gro cols
    setnames(
        grid_dt,
        c("seqnames", "start", "end", "i.width"),
        c("bin_chrom", "bin_start", "bin_end", "bin_width"),
        skip_absent = T
    )

    # adding ref and all to ideo backbone ranges not
    # present in the sv_dt
    grid_dt[is.na(cell), `:=`(
        cell = "all",
        sv_call_name = "ref",
        llr_to_ref = 0
    )]


    landscape_obj <- set_landscape_object(grid_dt, block)

    return(landscape_obj)
}
