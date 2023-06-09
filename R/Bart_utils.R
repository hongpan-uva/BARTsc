setMethod("str", "Bart", function(object) {
    cl <- class(object)
    # hide slots: 'Args', 'intermediate'
    sNms <- c("meta", "param", "data", "result") # slots to show
    # default arguments
    nest.lev <- 0
    give.head <- TRUE
    indent.str <- paste(rep.int(" ", max(0, nest.lev + 1)), collapse = "..")
    give.attr <- TRUE
    a <- attributes(object)

    trygetSlots <- function(x, nms) {
        r <- tryCatch(sapply(nms, methods::slot, object = x, simplify = FALSE),
            error = conditionMessage
        )
        if (is.list(r)) {
            r
        } else {
            warning("Not a validObject(): ", r, call. = FALSE) # instead of error
            r <- attributes(x) ## "FIXME" low-level assumption about S4 slots
            r <- r[names(r) != "class"]
            dp <- list(methods::getDataPart(x, NULL.for.none = TRUE))
            if (!is.null(dp)) names(dp) <- methods:::.dataSlot(nms)
            c(r, dp)
        }
    }

    strSub <- function(x, ...) {
        do.call(function(...) str(x, ...), c(list(), list(...)), quote = TRUE)
    }

    n.of. <- function(n, singl, plural) paste(n, ngettext(n, singl, plural))
    n.of <- function(n, noun) n.of.(n, noun, paste0(noun, "s"))

    cat("Formal class", " '", paste(cl, collapse = "', '"),
        "' [package \"", attr(cl, "package"), "\"] with ",
        n.of(length(sNms), "slot"), "\n",
        sep = ""
    )
    s <- trygetSlots(object, sNms)
    strSub(s,
        comp.str = "@ ", no.list = TRUE, give.length = give.head,
        indent.str = paste(indent.str, ".."), nest.lev = nest.lev + 1
    )
    return(invisible())
})


setGeneric("getResult", function(object, subcommand) standardGeneric("getResult"))

setMethod("getResult", "Bart", function(object, subcommand) {
    if (subcommand == "geneset") {
        result <- object@result$geneset$stats
    } else if (subcommand == "region") {
        result <- object@result$region$stats
    }

    result[, 2] <- round(result[, 2], 2)
    result[, 4] <- round(result[, 4], 3)
    result[, 5] <- round(result[, 5], 3)
    result[, 6] <- round(result[, 6], 3)

    return(result)
})
