#' combine the UDHS profile of modalities
#'
#' This function combine the UDHS profile of modalities to generate
#' a bimodal profile
#'
#' @param object a bart object
#'
#' @return a bart object
setGeneric("combineModalities", function(object) standardGeneric("combineModalities"))

setMethod("combineModalities", "bart", function(object) {
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

setMethod("getConsensusRank", "bart", function(object, method = "geom.mean") {
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

setMethod("combineModsByRank", "bart", function(object) {
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

#' combine the udhs profiles of two modalities based on rank aggregation
#' @param object a bart object
#'
#' @noRd
setGeneric("combineModsByTopRank", function(object, n_valid, method = "geom.mean") standardGeneric("combineModsByTopRank"))

setMethod("combineModsByTopRank", "bart", function(object, n_valid, method = "geom.mean") {
    if (!is.numeric(n_valid)) {
        warning("Input n_valid is not a numeric")
    }

    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])

    aggr_value <- .aggregate_rank(profile_1, profile_2, n_valid, method)

    # use the metric for aggregated ranking as 'counting'
    counting <- as.list(-aggr_value)

    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(object)
})

#' Aggregate two rank lists
#' 
#' @param vector_1 a numeric vector, names are entry ids
#' @param vector_2 a numeric vector, names are entry ids
#' @param n_valid number of top ranks to consider
#' @param method aggregation method
.aggregate_rank <- function(vector_1, vector_2, n_valid, method) {
    vector_1 <- vector_1[order(vector_1, decreasing = TRUE)]
    vector_2 <- vector_2[order(vector_2, decreasing = TRUE)]

    # set udhs ranked later than n_valid or non-positive to -Inf
    n_valid_pos_1 <- min(n_valid, length(which(vector_1 > 0))) + 1
    n_valid_pos_2 <- min(n_valid, length(which(vector_2 > 0))) + 1
    vector_1[n_valid_pos_1:length(vector_1)] <- -Inf
    vector_2[n_valid_pos_2:length(vector_2)] <- -Inf

    rank_1 <- rank(-vector_1, ties.method = "average")
    rank_2 <- rank(-vector_2, ties.method = "average")
    rank_1 <- rank_1[as.character(1:length(rank_1))]
    rank_2 <- rank_2[as.character(1:length(rank_2))]
    rank_df <- data.frame(rank_1, rank_2)

    if (method == "geom.mean") {
        geom_mean_vec <- apply(rank_df, 1, function(x) {
            return(exp(mean(log(x))))
        })
        return(geom_mean_vec)
    } else if (method == "MRR") {
        MRR_vec <- apply(rank_df, 1, function(x) {
            MRR <- ((1 / x[1]) + (1 / x[2])) / 2
            return(1 / MRR)
        })
        return(MRR_vec)
    } else if (method == "mean") {
        mean_vec <- apply(rank_df, 1, function(x) {
            return(mean(x))
        })
        return(mean_vec)
    }
}

#' combine the udhs profiles of two modalities based on the product of RNA and ATAC signals on UDHS
#' @param object a bart object
#' @param ATAC_weight weight of ATAC-seq
#'
#' @noRd
setGeneric("combineModsByProduct", function(object, ATAC_weight) standardGeneric("combineModsByProduct"))

setMethod("combineModsByProduct", "bart", function(object, ATAC_weight) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    id_num <- as.integer(names(profile_1))
    profile_1 <- profile_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    id_num <- as.integer(names(profile_2))
    profile_2 <- profile_2[order(id_num)]

    # profile_1[which(profile_1 < 0)] <- 0 # make RNA profile non-negative
    df <- data.frame(profile_1, profile_2)

    product <- df$profile_1 * (df$profile_2^ATAC_weight)
    # product <- (df$profile_1 + (df$profile_2^ATAC_weight)) / 2 # use mean profile
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

#' use the sign of RNA profile on ATAC profile
#' @param object a bart object
#'
#' @noRd
setGeneric("combineModsBySign", function(object) standardGeneric("combineModsBySign"))

setMethod("combineModsBySign", "bart", function(object) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    id_num <- as.integer(names(profile_1))
    profile_1 <- profile_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    id_num <- as.integer(names(profile_2))
    profile_2 <- profile_2[order(id_num)]

    df <- data.frame(profile_1, profile_2)

    df$profile_1[which(df$profile_1 > 0)] <- 1
    df$profile_1[which(df$profile_1 < 0)] <- -1

    profile_final <- df$profile_1 * df$profile_2
    names(profile_final) <- rownames(df)

    # use the profile_final as 'counting'
    counting <- as.list(profile_final)

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

setMethod("refineRNAWithATAC", "bart", function(object) {
    profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
    id_num <- as.integer(names(profile_1))
    profile_1 <- profile_1[order(id_num)]

    profile_2 <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
    id_num <- as.integer(names(profile_2))
    profile_2 <- profile_2[order(id_num)]

    # profile_1[which(profile_1 < 0)] <- 0 # make RNA profile non-negative
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
