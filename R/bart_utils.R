setMethod("str", "bart", function(object) {
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

#' return BART result
#'
#' This function return the result of BART
#'
#' @param object a bart object
#' @param subcommand "geneset" or "region" or "bimodal"
#'
#' @return a dataframe of bart result
setGeneric("get_bart_result", function(object, subcommand) standardGeneric("get_bart_result"))

setMethod("get_bart_result", "bart", function(object, subcommand) {
    if (subcommand == "geneset") {
        output.df <- object@result$geneset$stats
        output.df[, 2] <- round(output.df[, 2], 3)
        output.df[, 4] <- round(output.df[, 4], 3)
        output.df[, 5] <- round(output.df[, 5], 3)
        output.df[, 6] <- round(output.df[, 6], 3)
    } else if (subcommand == "region") {
        output.df <- object@result$region$stats
        output.df[, 2] <- round(output.df[, 2], 3)
        output.df[, 4] <- round(output.df[, 4], 3)
        output.df[, 5] <- round(output.df[, 5], 3)
        output.df[, 6] <- round(output.df[, 6], 3)
    } else if (subcommand == "bimodal") {
        output.df <- object@result$bimodal$stats
        output.df[, 2] <- round(output.df[, 2], 3)
        output.df[, 4] <- round(output.df[, 4], 3)
        output.df[, 5] <- round(output.df[, 5], 3)
        output.df[, 6] <- round(output.df[, 6], 3)
    }

    return(output.df)
})

#' return BART AUC values
#'
#' This function return the ROC AUC calculate by BART
#'
#' @param object a bart object
#' @param subcommand "geneset" or "region" or "bimodal"
#'
#' @return a dataframe of AUC values
setGeneric("get_bart_auc", function(object, subcommand) standardGeneric("get_bart_auc"))

setMethod("get_bart_auc", "bart", function(object, subcommand) {
    if (subcommand == "geneset") {
        output.df <- object@result$geneset$auc
    } else if (subcommand == "region") {
        output.df <- object@result$region$auc
    } else if (subcommand == "bimodal") {
        output.df <- object@result$bimodal$auc
    }

    return(output.df)
})


setGeneric("show_avail_result", function(object) standardGeneric("show_avail_result"))

setMethod("show_avail_result", "bart", function(object) {
    return(names(object@result))
})

setGeneric("show_avail_input", function(object) standardGeneric("show_avail_input"))

setMethod("show_avail_input", "bart", function(object) {
    return(names(object@data))
})
