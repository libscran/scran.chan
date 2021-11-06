transform_factor <- function(batch) {
    if (is.null(batch)) {
        list(index=NULL, names=NULL)
    } else {
        f <- factor(batch)
        list(index=as.integer(f) - 1L, names=levels(f))
    }
}

to_logical <- function(x, n, names) {
    if (is.logical(x)) {
        stopifnot(identical(length(x), n))
    } else if (is.numeric(x)) {
        stopifnot(min(x) >= 1 && max(x) <= n)
        tmp <- logical(n)
        tmp[x] <- TRUE
        x <- tmp
    } else if (is.character(x)) {
        stopifnot(!is.null(names))
        x <- names %in% x
    }
    x
}
