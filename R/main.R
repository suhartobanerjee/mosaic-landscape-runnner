library(stringr)

source("./class.R")
source("./set_grid.R")
source("./motion_calculations.R")
source("./penalise_on_sv.R")


get_sv <- function(object) {
    
    return(object@grid[object@current_bin, sv_call_name])
}


game_loop <- function(object) {

    if (object@current_bin < object@n_bins) {
        
        sv_bin <- get_sv(object)
        if(sv_bin != "ref") {
            if(str_detect(sv_bin, "dup")) {
                object <- penalise_dup(object)
            } else if(str_detect(sv_bin, "del")) {
                object <- penalise_del(object)
            } else if(str_detect(sv_bin, "inv")) {
                object <- penalise_inv(object)
            }
        }

        object@final_velocity <- calculate_final_velocity(object)
        object@time[length(object@time) + 1] <- calculate_time(object)

        # updating velocity and iteration counters
        object@initial_velocity <- object@final_velocity
        object@iter <- object@iter + 1
    }

    return(object)
}
#     time_vec <- cumsum(time_vec)

#     return(time_vec)
# }


# system.time(time <- game_loop(c(1:50000000)))
# time[1:10]
# cumsum(time)
# sum(time)
