am_i_ok <- function(ref, ptr, mutate=identity) {
    expect_identical(dim(ref), scran.chan:::tatami_dim(ptr$pointer))
    for (i in seq_len(ncol(ref))) {
        expected <- mutate(ref[,i])
        expect_identical(expected, scran.chan:::tatami_columns(ptr$pointer, i-1, 0, nrow(ref))[,1])
    }
}

compare <- function(ref, ptr) {
    am_i_ok(ref, ptr, mutate=unname)
}
