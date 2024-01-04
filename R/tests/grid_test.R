getwd()
setwd("../")
source("./main.R")

landscape_obj <- set_grid(block = 1e5)
landscape_obj

is_landscape_object(landscape_obj)
landscape_obj@grid[is.na(bin_chrom)]



# min(grid_dt$width)
# grid_dt[width == 0]
# grid_dt[id == 40]
