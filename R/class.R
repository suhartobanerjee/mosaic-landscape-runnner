setClass(
    Class = "Landscape",
    slots = c(
        grid_dt = "data.table",
        block = "numeric"
    )
)


is_landscape_object <- function(object) {
    # to check if landscape object

    if(class(object) == "Landscape") {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
