tatami_nrow <- function(x) {
    tatami_dim(x$pointer)[1]
}

tatami_ncol <- function(x) {
    tatami_dim(x$pointer)[2]
}

transform_factor <- function(f, n = NULL) {
    if (is.null(f)) {
        list(index=NULL, names=NULL)
    } else {
        if (!is.null(n)) {
            stopifnot(length(f) == n, !anyNA(f))
        }
        f <- factor(f)
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

do_downsample <- function(downsample) {
    !is.null(downsample) && downsample >= 1
}
