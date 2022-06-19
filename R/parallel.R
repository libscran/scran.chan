#' @importFrom parallel makeCluster
spawnCluster <- function(n) {
    if (.Platform$OS.type == "windows") {
        type <- "PSOCK"
    } else {
        type <- "FORK"
    }
    env <- new.env()
    env$cluster <- makeCluster(n, type=type)
    env$results <- list()
    env$active <- 0L
    env
}

.store_results <- function(env, value, tag) {
    type <- tag$type
    index <- tag$index

    if (! type %in% names(env$results)) {
        env$results[[type]] <- list()
    }

    # Appending it; we'll sort out the order later.    
    env$results[[type]] <- c(env$results[[type]], list(list(value = value, index = index)))
}

submitJob <- function(env, fun, args, type, index) {
    if (env$active == length(env$cluster)) {
        available <- parallel:::recvOneResult(env$cluster)
        .store_results(env, available$value, available$tag)
        to.use <- available$node
    } else {
        env$active <- env$active + 1L
        to.use <- env$active
    }

    tag <- list(type=type, index=index)
    parallel:::sendCall(env$cluster[[to.use]], fun, args, tag=tag)
}

finishJobs <- function(env) {
    for (i in seq_along(env$active)) {
        available <- parallel:::recvOneResult(env$cluster)
        .store_results(env, available$value, available$tag)
    }

    # Resorting everything.
    res <- env$results        
    for (type in names(res)) {
        current <- res[[type]]
        idx <- vapply(current, function(x) x$index, 0L)
        output <- vector("list", max(idx))
        output[idx] <- lapply(current, function(x) x$value)
        res[[type]] <- output
    }

    res
}
