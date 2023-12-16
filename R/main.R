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


game_loop <- function(rounds,
                      initial_velocity = 100,
                      acceleration = 50) {
   
    time_vec <- c(1:length(rounds))

    for (round in rounds) {
        final_velocity <- calculate_final_velocity(initial_velocity, acceleration)
        time_vec[round] <- calculate_time(initial_velocity, final_velocity)

        initial_velocity <- final_velocity
    }
    time_vec <- cumsum(time_vec)

    return(time_vec)
}


system.time(time <- game_loop(c(1:50000000)))
time[1:10]
cumsum(time)
sum(time)
