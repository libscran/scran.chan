#' @importFrom parallel makeCluster
spawnCluster <- function(n) {
    env <- new.env()

    # Don't attempt to be smart and use a FORK cluster on Unix-likes to squeeze
    # out more efficiency, as forking does not play nice with OpenMP; see
    # https://stackoverflow.com/questions/49049388/understanding-openmp-shortcomings-regarding-fork.
    # Just keep the PSOCK default.
    env$cluster <- makeCluster(n)

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
    tag <- list(type=type, index=index)

    if (is.null(env$cluster)) {
        res <- do.call(fun, args)
        .store_results(env, res, tag)

    } else {
        if (env$active == length(env$cluster)) {
            available <- parallel:::recvOneResult(env$cluster)
            .store_results(env, available$value, available$tag)
            to.use <- available$node
        } else {
            env$active <- env$active + 1L
            to.use <- env$active
        }
        parallel:::sendCall(env$cluster[[to.use]], fun, args, tag=tag)
    }
}

finishJobs <- function(env) {
    if (!is.null(env$cluster)) {
        for (i in seq_len(env$active)) {
            available <- parallel:::recvOneResult(env$cluster)
            .store_results(env, available$value, available$tag)
        }
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

#' @importFrom parallel stopCluster
.sweep_wrapper <- function(sweep, name, num.threads, drop) {
    params <- sweep(.env=NULL)

    if (nrow(params) == 1L || num.threads==1L) {
        env <- new.env()
        env$results <- list()
        env$cluster <- NULL
        env$active <- 0L
    } else {
        njobs <- nrow(params)
        nnodes <- min(njobs, num.threads)
        env <- spawnCluster(nnodes)
        num.threads <- max(1, floor(num.threads / nnodes))
        on.exit(stopCluster(env$cluster))
    }

    sweep(num.threads=num.threads, .env=env)
    completed <- finishJobs(env)

    if (nrow(params) == 1 && drop) {
        completed[[name]][[1]]
    } else {
        list(parameters = params, results = completed[[name]])
    }
}



