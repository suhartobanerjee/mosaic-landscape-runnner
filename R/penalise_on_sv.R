handle_dup <- function(object) {
    if (is_landscape_object(object)) {
        # first encounter of amp section
        if (object@amp_flag == 0) {
            object@amp_start_bin <- object@current_bin
            object@amp_flag <- object@amp_flag + 1
        }
    } else {
        stop("Landscape object not sent to handle function.")
    }

    object <- increment_current_bin(object)

    return(object)
}


handle_del <- function(object) {
    if (is_landscape_object(object)) {
        llr <- object@grids[[1]][
            bin_id == object@current_bin,
            llr_to_ref
        ]

        object@initial_velocity <- object@initial_velocity - object@initial_velocity / 2 * llr * object@prop_factor
    } else {
        stop("Landscape object not sent to handle function.")
    }

    object <- increment_current_bin(object)
    return(object)
}


handle_inv <- function(object) {
    if (is_landscape_object(object)) {
        llr <- object@grids[[1]][
            bin_id == object@current_bin,
            llr_to_ref
        ]

        object@acceleration <- object@acceleration - llr * object@prop_factor
    } else {
        stop("Landscape object not sent to handle function.")
    }

    object <- increment_current_bin(object)
    return(object)
}


handle_ref <- function(object) {
    if (str_detect(object@previous_sv, "dup") |
        is.na(object@previous_sv)) {
        # first repetition of the amp section
        if (object@amp_flag == 1) {
            object@current_bin <- object@amp_start_bin
            object@amp_start_bin <- NA_integer_
            object@amp_flag <- object@amp_flag + 1
        } else if (object@amp_flag == 2) {
            # resetting amp_flag
            object@amp_flag <- 0
        }
    }

    object <- increment_current_bin(object)
    return(object)
}
