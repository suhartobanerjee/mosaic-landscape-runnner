setClass(
    Class = "Landscape",
    slots = c(
        iter = "numeric",
        grid = "data.table",
        block = "numeric",
        n_bins = "numeric",
        initial_velocity = "numeric",
        final_velocity = "numeric",
        acceleration = "numeric",
        time = "numeric",
        prop_factor = "numeric",
        current_bin = "numeric",
        current_cell = "character",
        amp_start_bin = "numeric",
        amp_count = "numeric"
    )
)


set_landscape_object <- function(grid_dt, block) {
    
    landscape_obj <- new("Landscape",
        iter = 1,
        grid = grid_dt,
        block = block,
        n_bins = grid_dt[, max(bin_id)],
        initial_velocity = 100,
        final_velocity = NA_integer_,
        acceleration = 50,
        time = c(1),
        prop_factor = 0.1,
        current_bin = 1,
        current_cell = "",
        amp_start_bin = 0,
        amp_count = 0
    )

    return(landscape_obj)
}


is_landscape_object <- function(object) {
    # to check if landscape object

    if (class(object) == "Landscape") {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
