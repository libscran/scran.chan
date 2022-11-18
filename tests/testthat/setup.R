compare <- function(ref, ptr) {
    expect_identical(dim(ref), scran.chan:::tatami_dim(ptr$pointer))
    for (i in seq_len(ncol(ref))) {
        expect_identical(unname(ref[,i]), scran.chan:::tatami_columns(ptr$pointer, i-1, 0, nrow(ref))[,1])
    }
}
