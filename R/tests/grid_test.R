setwd("../")
source("./set_grid.R")

grid_dt <- set_grid()
grid_dt
min(grid_dt$width)
grid_dt[width == 0]
grid_dt[id == 40]
