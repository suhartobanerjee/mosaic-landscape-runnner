getwd()
setwd("../")
source("./main.R")

obj <- set_all_data(block = 1e5)
obj

# writing to file
fwrite(
    x = obj,
    file = "../data/tall3_proc.tsv.gz",
    sep = "\t"
)

all_cells <- obj@all_data[, unique(cell)]
all_cells

obj <- set_grid(obj, all_cells[2])
obj@cells
obj@grids
obj@grids[[1]][, unique(cell)]
obj


is_landscape_object(obj)
# landscape_obj@grid[is.na(bin_chrom)]



# min(grid_dt$width)
# grid_dt[width == 0]
# grid_dt[id == 40]
