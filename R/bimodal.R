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

#' get consensus rank of predictions results from both modalities
#' @param object a bart object
#' @param method aggregation method
#'
#' @noRd
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
#' @param list_1 the first ordered list to aggregate
#' @param list_2 the second ordered list to aggregate
#' @param method aggregation method
#'
#' @import RobustRankAggreg
#'
#' @noRd
.aggregate <- function(list_1, list_2, method = "geometric mean") {
    if (method %in% c("RRA", "min", "geom.mean", "mean", "median", "stuart")) {
        agg_result <- RobustRankAggreg::aggregateRanks(
            glist = list(list_1, list_2), method = method
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

#' combine the udhs profiles of two modalities based on rank aggregation
#' @param object a bart object
#'
#' @noRd
setGeneric("combineModsByRank", function(object) standardGeneric("combineModsByRank"))

setMethod("combineModsByRank", "Bart", function(object) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    rank_vector_1 <- rank(-profile_1, ties.method = "average")
    id_num <- as.integer(names(rank_vector_1))
    rank_vector_1 <- rank_vector_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    rank_vector_2 <- rank(-profile_2, ties.method = "average")
    id_num <- as.integer(names(rank_vector_2))
    rank_vector_2 <- rank_vector_2[order(id_num)]

    df <- data.frame(rank_vector_1, rank_vector_2)

    geom_mean_vec <- apply(df, 1, function(x) {
        return(exp(mean(log(x))))
    })

    # use the metric for aggregated ranking as 'counting'
    counting <- as.list(-geom_mean_vec)

    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(object)
})

#' combine the udhs profiles of two modalities based on the product of RNA and ATAC signals on UDHS
#' @param object a bart object
#'
#' @noRd
setGeneric("combineModsByProduct", function(object) standardGeneric("combineModsByProduct"))

setMethod("combineModsByProduct", "Bart", function(object) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    id_num <- as.integer(names(profile_1))
    profile_1 <- profile_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    id_num <- as.integer(names(profile_2))
    profile_2 <- profile_2[order(id_num)]

    df <- data.frame(profile_1, profile_2)

    product <- df$profile_1 * df$profile_2
    names(product) <- rownames(df)

    # use the product as 'counting'
    counting <- as.list(product)

    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(object)
})

#' combine the udhs profiles of two modalities based on the product of RNA and ATAC signals on UDHS
#' @param object a bart object
#'
#' @noRd
setGeneric("refineRNAWithATAC", function(object) standardGeneric("refineRNAWithATAC"))

setMethod("refineRNAWithATAC", "Bart", function(object) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    id_num <- as.integer(names(profile_1))
    profile_1 <- profile_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    id_num <- as.integer(names(profile_2))
    profile_2 <- profile_2[order(id_num)]

    df <- data.frame(profile_1, profile_2)

    df$profile_2 <- as.numeric(df$profile_2 > 0)

    product <- df$profile_1 * df$profile_2
    names(product) <- rownames(df)

    # use the product as 'counting'
    counting <- as.list(product)

    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(object)
})
