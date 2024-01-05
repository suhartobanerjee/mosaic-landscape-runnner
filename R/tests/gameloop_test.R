getwd()
setwd("../")
source("./main.R")
library(purrr)
library(ggplot2)


initialize_obj <- function(cell_id = 2) {
    obj <- set_all_data()
    obj <- set_grid(obj, obj@all_cells[cell_id])

    return(obj)
}
# obj@cells
# obj@grids
# obj@grids[[1]][, unique(cell)]
# obj@grids[[1]][, unique(sv_call_name)]
# obj
#
#
# obj@n_bins

obj <- initialize_obj(cell_id = 6)
obj

# ret_obj <- game_loop(obj)
# ret_obj
# ret_obj@iter

# obj@n_bins
# seq_along(1:obj@n_bins)
system.time(
    while (obj@current_bin < obj@n_bins) {
        obj <- game_loop(obj)
    }
)

obj@iter
nrow(obj@grids[[1]])
print(c(
    obj@initial_velocity,
    obj@final_velocity,
    obj@acceleration,
    obj@current_sv,
    obj@amp_start_bin,
    obj@amp_flag
))

# total time required
c_sum <- cumsum(obj@time)
c_sum[length(c_sum)]

# length(obj@time)

time_dt <- data.table(
    points = 1:length(obj@time),
    time = obj@time
)
time_dt

ggplot(
    data = time_dt,
    aes(
        x = points,
        y = time
    )
) +
    geom_line()

ggsave(str_glue("../plots/{obj@cells}_time_plots.pdf"))


obj@all_data[sv_call_name != "ref", .N, by = cell]
