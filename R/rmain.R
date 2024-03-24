library(stringr)

source("./class.R")
source("./set_grid.R")
source("./motion_calculations.R")
source("./penalise_on_sv.R")


increment_current_bin <- function(object) {
    prev_bin_idx <- object@grids[[1]][bin_id == object@current_bin, which = T]

    if (prev_bin_idx < nrow(object@grids[[1]])) {
        object@current_bin <- object@grids[[1]][prev_bin_idx + 1, bin_id]
    }

    return(object)
}


set_prev_curr_sv <- function(object) {
    object@previous_sv <- object@current_sv
    object@current_sv <- object@grids[[1]][
        bin_id == object@current_bin,
        sv_call_name
    ]

    return(object)
}


game_loop <- function(object) {
    object <- set_prev_curr_sv(object)

    # handle gen status and increment current_bin
    if (object@current_sv != "ref") {
        if (str_detect(object@current_sv, "dup")) {
            object <- handle_dup(object)
        } else if (str_detect(object@current_sv, "del")) {
            object <- handle_del(object)
        } else if (str_detect(object@current_sv, "inv")) {
            object <- handle_inv(object)
        }
    } else {
        object <- handle_ref(object)
    }

    # motion_calculations
    object@final_velocity <- calculate_final_velocity(object)
    object@time[length(object@time) + 1] <- calculate_time(object)

    # updating velocity and iteration counters
    object@initial_velocity <- object@final_velocity
    object@iter <- object@iter + 1


    return(object)
}
