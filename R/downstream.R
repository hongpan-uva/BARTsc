setGeneric("find_differential_TFs", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    standardGeneric("find_differential_TFs")
})

setMethod("find_differential_TFs", "bartsc", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    if (min_diff > 1 || min_diff < -1) {
        stop("min_diff is supposed to be between -1 and 1")
    }

    if (missing(mod)) {
        mod <- object@meta$active_mod
    }

    if (missing(cell_types_used)) {
        cell_types_used <- object@meta$cell_types_used
    }

    matrix_list <- object@resultsCrossCellType[[mod]]$deviation

    cds_df <- object@resultsCrossCellType[[mod]]$CDS
    tf_used <- cds_df$TF[which(cds_df$n_profiles >= min_N_profile)]

    pos_tf_list <- c()
    neg_tf_list <- c()

    if (missing(group2)) {
        group2 <- cell_types_used[which(!cell_types_used %in% group1)]
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (all(mtx[group1, group2] >= min_diff)) {
            pos_tf_list <- c(pos_tf_list, tf)
        }
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (all(mtx[group1, group2] <= -min_diff)) {
            neg_tf_list <- c(neg_tf_list, tf)
        }
    }

    return(list(positive = pos_tf_list, negative = neg_tf_list))
})


setGeneric("find_differential_TFs_2", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    standardGeneric("find_differential_TFs_2")
})

setMethod("find_differential_TFs_2", "bartsc", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    if (min_diff > 1 || min_diff < -1) {
        stop("min_diff is supposed to be between -1 and 1")
    }

    if (missing(mod)) {
        mod <- object@meta$active_mod
    }

    if (missing(cell_types_used)) {
        cell_types_used <- object@meta$cell_types_used
    }

    matrix_list <- object@resultsCrossCellType[[mod]]$deviation

    cds_df <- object@resultsCrossCellType[[mod]]$CDS
    tf_used <- cds_df$TF[which(cds_df$n_profiles >= min_N_profile)]

    pos_tf_list <- c()
    neg_tf_list <- c()

    if (missing(group2)) {
        group2 <- cell_types_used[which(!cell_types_used %in% group1)]
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (mean(mtx[group1, group2]) >= min_diff) {
            pos_tf_list <- c(pos_tf_list, tf)
        }
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (mean(mtx[group1, group2]) <= -min_diff) {
            neg_tf_list <- c(neg_tf_list, tf)
        }
    }

    return(list(positive = pos_tf_list, negative = neg_tf_list))
})

setGeneric("find_differential_TFs_3", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    standardGeneric("find_differential_TFs_3")
})


setMethod("find_differential_TFs_3", "bartsc", function(
    object,
    cell_types_used,
    group1,
    group2,
    mod,
    min_diff = 0,
    min_N_profile = 3) {
    if (min_diff > 1 || min_diff < -1) {
        stop("min_diff is supposed to be between -1 and 1")
    }

    if (missing(mod)) {
        mod <- object@meta$active_mod
    }

    if (missing(cell_types_used)) {
        cell_types_used <- object@meta$cell_types_used
    }

    matrix_list <- object@resultsCrossCellType[[mod]]$deviation

    cds_df <- object@resultsCrossCellType[[mod]]$CDS
    tf_used <- cds_df$TF[which(cds_df$n_profiles >= min_N_profile)]

    pos_tf_list <- c()
    neg_tf_list <- c()

    if (missing(group2)) {
        group2 <- cell_types_used[which(!cell_types_used %in% group1)]
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (mean(mtx[group1, group2]) >= min_diff) {
            pos_tf_list <- c(pos_tf_list, tf)
        }
    }

    for (tf in tf_used) {
        mtx <- matrix_list[[tf]]
        if (mean(mtx[group1, group2]) <= -min_diff) {
            neg_tf_list <- c(neg_tf_list, tf)
        }
    }

    return(list(positive = pos_tf_list, negative = neg_tf_list))
})

#' Make a 3d scatter plot for key regulators identification
#'
#' BARTsc identifies key regulators for a given cell types by comprehensively
#' considering three factors: 1) Contribution to signature expression, 2) rela-
#' tive activity across all the cell types and 3) the uniqueness of a given
#' transcription regulator. These three factors are quantified by signature
#' score, MDR (mean deviation ratio) and dMDR (differential MDR) respectively.
#' Final rank is determined by the mean rank of above indices. Final pvalues
#' are calculated based on Irvin-Hall distribution. To garantee the robustness
#' of final output, only transcription regulators with at least 3 binding pro-
#' files in database are considered by default.
#'
#' @param object A bartsc object.
#' @param mod One of "RNA", "ATAC" and "bimodal".
#' @param min.N.profile Only transcription regulators with at least certain
#' number of binding profiles in database are considered. Default is 3.
#'
#' @return A bartsc object.
#'
#' @export
#'
setGeneric("find_key_regulators", function(
    object,
    mod,
    min.N.profile = 3) {
    standardGeneric("find_key_regulators")
})

setMethod("find_key_regulators", "bartsc", function(
    object,
    mod,
    min.N.profile = 3) {
    matrices <- object@resultsCrossCellType[[mod]]$deviation

    # calculate mean DR
    mdr_list <- lapply(matrices, function(x) {
        diag(x) <- NA
        return(rowMeans(x, na.rm = TRUE))
    })

    # filter tfs by profile numbers
    cds_df <- object@resultsCrossCellType[[mod]]$CDS
    tfs_used <- cds_df[which(cds_df$n_profiles >= min.N.profile), "TF"]

    # rank aggregation by mean
    RankAggr.mean <- function(data, cols, weight, return.rank = TRUE) {
        ncol <- length(cols)
        if (missing(weight)) {
            weight <- rep(1, ncol)
        }

        rank_vec <- c()
        for (i in 1:ncol) {
            cl <- cols[i]
            rank_vec <- c(rank_vec, rank(-data[, cl], ties.method = "average") * weight[i])
        }

        rank.df <- as.data.frame(matrix(rank_vec, ncol = length(cols)))

        rank.df$aggr.value <- rowMeans(rank.df)

        rank.df$final_rank <- rank(rank.df$aggr.value)
        if (return.rank) {
            return(rank.df$final_rank)
        } else {
            return(rank.df)
        }
    }

    # irwinhall distribution
    cdf_irwinhall <- function(x, n) {
        if (x < 0 | x > n) {
            stop("x should be within [0,n]")
        }

        sum <- 0
        k <- 0
        while (k <= floor(x)) {
            sum <- sum + (-1)**k * choose(n, k) * (x - k)**n
            k <- k + 1
        }

        cdf <- sum / factorial(n)
        return(cdf)
    }

    output_list <- list()

    for (ct in object@meta$cell_types_used) {
        sig_res <- get_result(object, analysis = "cell type signature")[[ct]]
        sig_res$signature_score <- -log10(sig_res$rank_avg_z_p_a_irwinhall_pvalue)
        output_df <- sig_res[, c(1, 9)]

        # MDR
        MDR_vec <- c()
        for (i in output_df$TF) {
            MDR_vec <- c(MDR_vec, mdr_list[[i]][ct])
        }
        output_df$MDR <- MDR_vec

        # difference of MDR between the target cell type and the next
        # highest one
        dMDR_vec <- c()
        for (i in output_df$TF) {
            tmp <- mdr_list[[i]]
            v1 <- tmp[ct]
            v2 <- max(tmp[which(names(tmp) != ct)])
            dMDR_vec <- c(dMDR_vec, v1 - v2)
        }
        output_df$dMDR <- dMDR_vec

        # filter by tf used
        output_df <- output_df[which(output_df$TF %in% tfs_used), ]

        rank_res <- RankAggr.mean(output_df, c("signature_score", "MDR", "dMDR"), c(1, 1, 1), return.rank = FALSE)
        output_df <- cbind(output_df, rank_res[, c("aggr.value", "final_rank")])

        # calculate pvalue
        output_df$final_pvalue <- sapply(output_df$aggr.value, function(x) {
            return(cdf_irwinhall(3 * x / nrow(output_df), 3))
        })
        output_df <- output_df[order(output_df$final_rank, decreasing = F), ]
        rownames(output_df) <- 1:nrow(output_df)

        output_list[[ct]] <- output_df
    }

    object@resultsKeyRegsIdent[[mod]] <- output_list
    return(object)
})
