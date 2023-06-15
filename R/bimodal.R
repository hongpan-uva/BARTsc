#' combine the UDHS profile of modalities
#'
#' This function combine the UDHS profile of modalities to generate
#' a bimodal profile
#'
#' @param object a Bart object
#'
#' @return a Bart object
setGeneric("combineModalities", function(object) standardGeneric("combineModalities"))

setMethod("combineModalities", "Bart", function(object) {
    gene_score_li <- object@intermediate[["Marge_based"]][["predicted_enhancers"]]
    gene_score_df <- data.frame(
        ID = as.integer(names(gene_score_li)),
        gene_score = as.numeric(unname(unlist(gene_score_li)))
    ) %>%
        dplyr::arrange(., ID)

    region_score_li <- object@intermediate[["region_based"]][["overlapped_enhancers"]]
    region_score_df <- data.frame(
        ID = as.integer(names(region_score_li)),
        region_score = as.numeric(unname(unlist(region_score_li)))
    ) %>%
        dplyr::arrange(., ID)

    bimodal_score <- cbind(
        gene_score_df[, c("ID", "gene_score")],
        region_score_df[, c("region_score")]
    )
    colnames(bimodal_score) <- c("ID", "gene_score", "region_score")

    # combine scores
    bimodal_score$combined_score <- apply(bimodal_score, 1, function(x) {
        return(sum(x[c(2, 3)]))
    })

    counting <-
        bimodal_score$combined_score %>%
        split(., bimodal_score$ID)

    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(object)
})


setGeneric("getConsensusRank", function(object, method = "geom.mean") standardGeneric("getConsensusRank"))

setMethod("getConsensusRank", "Bart", function(object, method = "geom.mean") {
    ava_results <- showAvailableResult(object)
    if (!"geneset" %in% ava_results) {
        stop("BART result of geneset is missing")
    } else if (!"region" %in% ava_results) {
        stop("BART result of region is missing")
    }

    gene_res <- getResult(object, "geneset")
    region_res <- getResult(object, "region")

    if (!setequal(gene_res$TF, region_res$TF)) {
        stop("modalities to aggregate have different TFs")
    }

    rank_table <- .aggregate(gene_res$TF, region_res$TF, method = method)

    if ("bimodal" %in% ava_results) { # remove the existing result
        object@result[["bimodal"]] <- NULL
    }
    object@result[["bimodal"]][["rank_table"]] <- rank_table

    return(object)
})

#' aggregate rank of TFs from different modalities
#' @param rank_df a dataframe, 2nd and 3rd col: first and second TF ranks
#' @param method aggregation method
#'
#' @import RobustRankAggreg
#'
#' @noRd
.aggregate <- function(list_1, list_2, method = "geometric mean") {
    if (method %in% c("RRA", "min", "geom.mean", "mean", "median", "stuart")) {
        agg_result <- RobustRankAggreg::aggregateRanks(
            glist = list(list_1, list_2), method = "median"
        )

        rank_1 <- sapply(agg_result$Name, function(x) {
            which(x == list_1)
        })
        rank_2 <- sapply(agg_result$Name, function(x) {
            which(x == list_2)
        })

        output.df <- data.frame(
            TF = agg_result$Name,
            rank_geneset = rank_1,
            rank_region = rank_2,
            consensus_rank = seq_len(nrow(agg_result)),
            score = agg_result$Score
        )

        colnames(output.df)[5] <- paste0(method, "_score")
        rownames(output.df) <- agg_result$Name
    }

    return(output.df)
}
