calculate_final_velocity <- function(object) {
    final_velocity <- sqrt(object@initial_velocity^2 + 2 * object@acceleration * object@block)

    return(final_velocity)
}


calculate_time <- function(object) {
    time <- 2 * object@block / (object@initial_velocity + object@final_velocity)

    return(time)
}
