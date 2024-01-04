penalise_dup <- function(object) {

    if (is_landscape_object(object)) {

        # first encounter of amp section
        if(object@amp_count == 0) {
            object@amp_start_bin <- object@current_bin
            objec@amp_count <- object@amp_count + 1

            return(object)
        }


    } else {
        stop("Landscape object not sent to penalise function.")
    }

    return(object)
}


penalise_del <- function(object) {

    if (is_landscape_object(object)) {
        prop_factor <- object@grid[
            cell == object@current_cell &
                bin_id == object@current_bin,
            llr_to_ref
        ]

        object@acceleration <- object@acceleration - object@acceleration * prop_factor
    } else {
        stop("Landscape object not sent to penalise function.")
    }

    return(object)
}


penalise_inv <- function(object) {

    if (is_landscape_object(object)) {
        prop_factor <- object@grid[
            cell == object@current_cell &
                bin_id == object@current_bin,
            llr_to_ref
        ]

        object@acceleration <- object@acceleration - object@acceleration * prop_factor
    } else {
        stop("Landscape object not sent to penalise function.")
    }

    return(object)
}
