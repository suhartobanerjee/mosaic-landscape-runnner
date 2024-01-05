setClass(
    Class = "Landscape",
    slots = c(
        iter = "numeric",
        all_data = "data.table",
        all_cells = "character",
        grids = "list",
        block = "numeric",
        n_bins = "numeric",
        initial_velocity = "numeric",
        final_velocity = "numeric",
        acceleration = "numeric",
        time = "numeric",
        prop_factor = "numeric",
        current_bin = "numeric",
        cells = "character",
        previous_sv = "character",
        current_sv = "character",
        amp_start_bin = "numeric",
        amp_flag = "numeric"
    )
)


set_landscape_object <- function(all_data, block) {
    landscape_obj <- new("Landscape",
        iter = 1,
        all_data = all_data,
        all_cells = all_data[, unique(cell)],
        grids = list(),
        block = block,
        n_bins = NA_integer_,
        initial_velocity = 100,
        final_velocity = NA_integer_,
        acceleration = 50,
        time = numeric(0),
        prop_factor = 0.01,
        current_bin = NA_integer_,
        cells = character(0),
        previous_sv = NA_character_,
        current_sv = NA_character_,
        amp_start_bin = 0,
        amp_flag = 0
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
