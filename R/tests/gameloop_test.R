getwd()
setwd("../")
source("./main.R")
library(purrr)
library(ggplot2)


obj <- set_grid()
obj

obj@grid <- obj@grid[cell == "TALL3x1_DEA5_PE20419" |
         cell == "all"]

obj@n_bins
obj@grid[sv_call_name == "dup_h2"]
# ret_obj <- game_loop(obj)
# ret_obj


for(x in c(1:obj@n_bins)) {
    obj@current_bin <- x
    obj <- game_loop(obj)
}
obj@iter
obj@initial_velocity
obj@final_velocity

obj@amp_start_bin
obj@amp_count

obj@time
length(obj@time)
obj@time[167]

time_dt <- data.table(
    points = 1:obj@n_bins,
    time = obj@time[1:obj@n_bins]
)
time_dt

ggplot(data = time_dt,
       aes(x = points,
           y = time
       )
) +
    geom_line()

ggsave("../plots/time_plots.pdf")
