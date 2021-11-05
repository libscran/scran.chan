transform_factor <- function(batch) {
    if (is.null(batch)) {
        list(index=NULL, names=NULL)
    } else {
        f <- factor(batch)
        list(index=as.integer(f) - 1L, names=levels(f))
    }
}
