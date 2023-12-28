calculate_final_velocity <- function(initial_velocity, acceleration) {
   
    distance <- 200
    final_velocity <- sqrt(initial_velocity ^ 2 + 2 * acceleration * distance)

    return(final_velocity)
}


calculate_time <- function(initial_velocity, final_velocity) {
    
    distance <- 200
    time <- 2 * distance / (initial_velocity + final_velocity)

    return(time)
}


