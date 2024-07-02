#' @importFrom Matrix t
calc_pct_fc_rna <- function(matrix, group1, group2) {
    # calculate percentage and foldchange

    listCols <- function(m) {
        # converts a sparse Matrix into a list of its columns
        res <- split(m@x, findInterval(seq_len(length(m@x)), m@p, left.open = TRUE))

        # get numbers of non-zero elements of each columns
        col_nnzero <- m@p[2:length(m@p)] - m@p[1:(length(m@p) - 1)]
        names(res) <- colnames(m)[(which(col_nnzero > 0))]

        # for the columns with all 0, create an empty list fot them
        to_append <- sapply(which(col_nnzero == 0), function(x) {
            return(list())
        })
        names(to_append) <- colnames(m)[(which(col_nnzero == 0))]

        # append two lists and re-order by colname
        res_final <- append(res, to_append)
        res_final <- res_final[colnames(m)]

        return(res_final)
    }

    start <- Sys.time()
    fore_mtx <- t(matrix[, group1])
    back_mtx <- t(matrix[, group2])

    list_fore <- listCols(fore_mtx)
    list_back <- listCols(back_mtx)

    res_df <- sapply(1:nrow(matrix), function(x) {
        fore_vec <- list_fore[[x]]
        back_vec <- list_back[[x]]
        fore_pct <- length(fore_vec) / length(group1)
        back_pct <- length(back_vec) / length(group2)
        if (length(fore_vec) == 0) {
            data.1 <- 0
        } else {
            data.1 <- log2((sum(expm1(fore_vec)) / length(group1)) + 1)
        }
        if (length(back_vec) == 0) {
            data.2 <- 0
        } else {
            data.2 <- log2((sum(expm1(back_vec)) / length(group2)) + 1)
        }

        log2fc <- data.1 - data.2

        return(c(log2fc, fore_pct, back_pct, data.1, data.2))
    })

    end <- Sys.time()
    tcost <- end - start
    print(paste("Pct and FC calculation:", round(as.numeric(tcost), 4), units(tcost)))

    res_df <- as.data.frame(t(res_df))
    colnames(res_df) <- c("log2fc", "fore_pct", "back_pct", "fore_avg_exp", "back_avg_exp")
    rownames(res_df) <- rownames(matrix)

    return(res_df)
}

calc_pct_fc_atac <- function(matrix, group1, group2) {
    # calculate percentage and foldchange

    listCols <- function(m) {
        # converts a sparse Matrix into a list of its columns
        res <- split(m@x, findInterval(seq_len(length(m@x)), m@p, left.open = TRUE))

        # get numbers of non-zero elements of each columns
        col_nnzero <- m@p[2:length(m@p)] - m@p[1:(length(m@p) - 1)]
        names(res) <- colnames(m)[(which(col_nnzero > 0))]

        # for the columns with all 0, create an empty list fot them
        to_append <- sapply(which(col_nnzero == 0), function(x) {
            return(list())
        })
        names(to_append) <- colnames(m)[(which(col_nnzero == 0))]

        # append two lists and re-order by colname
        res_final <- append(res, to_append)
        res_final <- res_final[colnames(m)]

        return(res_final)
    }

    start <- Sys.time()
    fore_mtx <- t(matrix[, group1])
    back_mtx <- t(matrix[, group2])

    list_fore <- listCols(fore_mtx)
    list_back <- listCols(back_mtx)

    res_df <- sapply(1:nrow(matrix), function(x) {
        fore_vec <- list_fore[[x]]
        back_vec <- list_back[[x]]
        fore_pct <- length(fore_vec) / length(group1)
        back_pct <- length(back_vec) / length(group2)
        if (length(fore_vec) == 0) {
            sum_fore <- 0
        } else {
            sum_fore <- sum(fore_vec)
        }
        if (length(back_vec) == 0) {
            sum_back <- 0
        } else {
            sum_back <- sum(back_vec)
        }
        fore_avg_exp <- sum_fore / length(group1)
        back_avg_exp <- sum_back / length(group2)
        log2fc <- log2((fore_avg_exp + 1) / (back_avg_exp + 1))
        # log2fc <- log2((sum_fore * length(group2) + 1) / (sum_back * length(group1) + 1))
        return(c(log2fc, fore_pct, back_pct, fore_avg_exp, back_avg_exp))
    })

    end <- Sys.time()
    tcost <- end - start
    print(paste("Pct and FC calculation:", round(as.numeric(tcost), 4), units(tcost)))

    res_df <- as.data.frame(t(res_df))
    colnames(res_df) <- c("log2fc", "fore_pct", "back_pct", "fore_avg_exp", "back_avg_exp")

    return(res_df)
}

#' find differentially expressed features
#'
#' @param object A bartsc object.
#'
#' @importFrom presto wilcoxauc
#'
#' @return A character vector.
setGeneric("find_pos_diffs", function(
    object,
    mod,
    fore_celltype,
    back_celltype,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    random.seed = 1) {
    standardGeneric("find_pos_diffs")
})

setMethod("find_pos_diffs", "bartsc", function(
    object,
    mod,
    fore_celltype,
    back_celltype,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    random.seed = 1) {
    if (missing(mod)) {
        stop("must specify a valid modality ('RNA' or 'ATAC')")
    } else if (!mod %in% c("RNA", "ATAC")) {
        stop("must specify a valid modality ('RNA' or 'ATAC')")
    }
    if (mod == "RNA") {
        mtx_use <- object@assays$RNA$normalized
    } else if (mod == "ATAC") {
        mtx_use <- object@assays$ATAC$normalized
    }

    label <- object@meta$label
    fore_idx <- which(label %in% fore_celltype)
    back_idx <- which(label %in% back_celltype)

    # calculate percentage and foldchange
    if (mod == "RNA") {
        sum_df <- calc_pct_fc_rna(mtx_use, fore_idx, back_idx)
    } else if (mod == "ATAC") {
        sum_df <- calc_pct_fc_atac(mtx_use, fore_idx, back_idx)
    }

    # label2 for presto wilcoxon test
    label2 <- object@meta$label
    levels(label2)[which(levels(label2) %in% fore_celltype)] <- "fore"
    levels(label2)[which(levels(label2) %in% back_celltype)] <- "back"

    # presto wilcoxon test
    wilcox_df <- presto::wilcoxauc(mtx_use, groups_use = c("fore", "back"), label2)
    wilcox_df <- wilcox_df[which(wilcox_df$group == "fore"), ]

    final_de_df <- cbind(sum_df, wilcox_df[, c("feature", "auc", "pval", "padj")])

    return(final_de_df)
})

#' Find cell type signature genes
#'
#' Find cell type signature genes by comparing the target cell type with other
#' cells. Wilcoxon rank sum test and Benjamini-Hochberg FDR adjustment.
#'
#' @param object A bartsc object.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the
#' function by not testing genes that are very infrequently expressed. Default
#' is 0.1.
#' @param min.diff.pct Only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default.
#' @param log2fc.thr Threshold of log2(foldchange).
#' @param pval.thr Threshold of p value.
#' @param padj.thr Threshold of adjusted p value.
#' @param auc.thr Threshold of AUC calculated by presto.
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf).
#'
#' @return object A bartsc object.
#'
#' @export
#'
setGeneric("find_celltype_deg", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    standardGeneric("find_celltype_deg")
})

#' @import parallel
setMethod("find_celltype_deg", "bartsc", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    object@data[["signature_genes"]] <- list()
    cell_types_used <- object@meta$cell_types_used
    object@data[["signature_genes"]] <- lapply(
        X = cell_types_used,
        FUN = function(x) {
            back_celltype <- cell_types_used[!cell_types_used == x]
            DE_df <- find_pos_diffs(
                object = object,
                mod = "RNA",
                fore_celltype = x,
                back_celltype = back_celltype,
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                max.cells.per.ident = max.cells.per.ident
            )
            if (!is.null(pval.thr)) {
                DE_df <- DE_df[which(DE_df$pval < pval.thr), ]
            }
            if (!is.null(padj.thr)) {
                DE_df <- DE_df[which(DE_df$padj < padj.thr), ]
            }
            DE_df <- DE_df[
                DE_df$log2fc >= log2fc.thr &
                    DE_df$fore_pct >= min.pct &
                    DE_df$auc >= auc.thr,
            ]
            message(paste0(x, ":", nrow(DE_df)))
            return(DE_df$feature)
        }
    )
    names(object@data[["signature_genes"]]) <- cell_types_used
    return(object)
})

#' Find cell type signature peaks
#'
#' Find cell type signature peaks by comparing the target cell type with other
#' cells. Wilcoxon rank sum test and Benjamini-Hochberg FDR adjustment.
#'
#' @param object A bartsc object.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the
#' function by not testing genes that are very infrequently expressed. Default
#' is 0.1.
#' @param min.diff.pct Only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default.
#' @param log2fc.thr Threshold of log2(foldchange).
#' @param pval.thr Threshold of p value.
#' @param padj.thr Threshold of adjusted p value.
#' @param auc.thr Threshold of AUC calculated by presto.
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf).
#'
#' @return object A bartsc object.
#'
#' @export
#'
setGeneric("find_celltype_dar", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    standardGeneric("find_celltype_dar")
})

setMethod("find_celltype_dar", "bartsc", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    object@data[["signature_peaks"]] <- list()
    cell_types_used <- object@meta$cell_types_used

    object@data[["signature_peaks"]] <- lapply(
        X = cell_types_used,
        FUN = function(x) {
            back_celltype <- cell_types_used[!cell_types_used == x]
            DE_df <- find_pos_diffs(
                object = object,
                mod = "ATAC",
                fore_celltype = x,
                back_celltype = back_celltype,
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                max.cells.per.ident = max.cells.per.ident
            )
            if (!is.null(pval.thr)) {
                DE_df <- DE_df[which(DE_df$pval < pval.thr), ]
            }
            if (!is.null(padj.thr)) {
                DE_df <- DE_df[which(DE_df$padj < padj.thr), ]
            }
            DE_df <- DE_df[
                DE_df$log2fc >= log2fc.thr &
                    DE_df$fore_pct >= min.pct &
                    DE_df$auc >= auc.thr,
            ]
            out_df <- object@assays$ATAC$peaks[DE_df$feature, ]
            out_df$score <- DE_df$log2fc
            message(paste0(x, ":", nrow(out_df)))
            return(out_df)
        }
    )

    names(object@data[["signature_peaks"]]) <- cell_types_used
    return(object)
})

setGeneric("find_celltype_ar", function(object, pct = 0.05) {
    standardGeneric("find_celltype_ar")
})

setMethod("find_celltype_ar", "bartsc", function(object, pct = 0.05) {
    if (pct < 0 || pct > 1) {
        stop("invalid pct argument")
    }

    message("Calling cell type accessible regions...")

    norm_mtx <- object@assays$ATAC$normalized
    options(scipen = 200)

    object@data$accessible_peaks <- list()

    cell_types_used <- object@meta$cell_types_used
    cell_type_label <- object@meta$label
    for (ct in cell_types_used) {
        cells <- names(cell_type_label)[which(cell_type_label == ct)]
        threshold <- length(cells) * pct
        acc_score <- apply(norm_mtx[, cells], 1, function(x) {
            if (length(x[which(x > 0)]) > threshold) {
                return(mean(x))
            } else {
                return(0)
            }
        })
        all_peaks <- object@assays$ATAC$peaks
        all_peaks$score <- acc_score
        acc_peaks <- all_peaks[which(all_peaks$score > 0), ]
        object@data$accessible_peaks[[ct]] <- acc_peaks
        message(paste(ct, nrow(acc_peaks)))
    }

    return(object)
})

#' Find pairwise differentially expressed genes
#'
#' Find cell type differentially expressed genes between each pair of cell
#' types. Wilcoxon rank sum test and Benjamini-Hochberg FDR adjustment.
#'
#' @param object A bartsc object.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the
#' function by not testing genes that are very infrequently expressed. Default
#' is 0.1.
#' @param min.diff.pct Only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default.
#' @param log2fc.thr Threshold of log2(foldchange).
#' @param pval.thr Threshold of p value.
#' @param padj.thr Threshold of adjusted p value.
#' @param auc.thr Threshold of AUC calculated by presto.
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf).
#'
#' @return object A bartsc object.
#'
#' @export
#'
setGeneric("find_pairwise_deg", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    standardGeneric("find_pairwise_deg")
})

setMethod("find_pairwise_deg", "bartsc", function(
    object,
    min.pct = 0.1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    cell_types_used <- object@meta$cell_types_used
    pairs <- list()
    pair_names <- c()
    idx <- 1

    for (ct1 in cell_types_used) {
        for (ct2 in cell_types_used) {
            if (ct1 == ct2) {
                next()
            }
            pairs[[idx]] <- c(ct1, ct2)
            pair_names <- c(pair_names, paste0(ct1, "::", ct2))
            idx <- idx + 1
        }
    }

    message("Calling cell-type pairwise DEG...")

    object@data[["pairwise_DEG"]] <- lapply(
        X = pairs,
        FUN = function(x) {
            DE_df <- find_pos_diffs(
                object = object,
                mod = "RNA",
                fore_celltype = x[1],
                back_celltype = x[2],
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                max.cells.per.ident = max.cells.per.ident
            )
            if (!is.null(pval.thr)) {
                DE_df <- DE_df[which(DE_df$pval < pval.thr), ]
            }
            if (!is.null(padj.thr)) {
                DE_df <- DE_df[which(DE_df$padj < padj.thr), ]
            }
            DE_df <- DE_df[
                DE_df$log2fc >= log2fc.thr &
                    DE_df$fore_pct >= min.pct &
                    DE_df$auc >= auc.thr,
            ]
            message(paste0(x[1], "::", x[2], ": ", nrow(DE_df)))
            return(DE_df$feature)
        }
    )

    names(object@data[["pairwise_DEG"]]) <- pair_names
    return(object)
})

#' Find pairwise differentially accessible regions
#'
#' Find cell type differentially accessible regions between each pair of cell
#' types. Wilcoxon rank sum test and Benjamini-Hochberg FDR adjustment.
#'
#' @param object A bartsc object.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the
#' function by not testing genes that are very infrequently expressed. Default
#' is 0.1.
#' @param min.diff.pct Only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default.
#' @param log2fc.thr Threshold of log2(foldchange).
#' @param pval.thr Threshold of p value.
#' @param padj.thr Threshold of adjusted p value.
#' @param auc.thr Threshold of AUC calculated by presto.
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf).
#'
#' @return object A bartsc object.
#'
#' @export
#'
setGeneric("find_pairwise_dar", function(
    object,
    min.pct = 0.1,
    max.back.pct = 1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    standardGeneric("find_pairwise_dar")
})

setMethod("find_pairwise_dar", "bartsc", function(
    object,
    min.pct = 0.1,
    max.back.pct = 1,
    min.diff.pct = -Inf,
    log2fc.thr = 0.25,
    pval.thr = NULL,
    padj.thr = NULL,
    auc.thr = 0,
    max.cells.per.ident = Inf) {
    cell_types_used <- object@meta$cell_types_used
    pairs <- list()
    pair_names <- c()
    idx <- 1

    for (ct1 in cell_types_used) {
        for (ct2 in cell_types_used) {
            if (ct1 == ct2) {
                next()
            }
            pairs[[idx]] <- c(ct1, ct2)
            pair_names <- c(pair_names, paste0(ct1, "::", ct2))
            idx <- idx + 1
        }
    }

    message("Calling cell-type pairwise DAR...")

    object@data[["pairwise_DAR"]] <- lapply(
        X = pairs,
        FUN = function(x) {
            DE_df <- find_pos_diffs(
                object = object,
                mod = "ATAC",
                fore_celltype = x[1],
                back_celltype = x[2],
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                max.cells.per.ident = max.cells.per.ident
            )
            if (!is.null(pval.thr)) {
                DE_df <- DE_df[which(DE_df$pval < pval.thr), ]
            }
            if (!is.null(padj.thr)) {
                DE_df <- DE_df[which(DE_df$padj < padj.thr), ]
            }
            DE_df <- DE_df[
                DE_df$log2fc >= log2fc.thr &
                    DE_df$fore_pct >= min.pct &
                    DE_df$auc >= auc.thr,
            ]
            message(paste0(x[1], "::", x[2], ": ", nrow(DE_df)))
            out_df <- object@assays$ATAC$peaks[DE_df$feature, ]
            out_df$score <- DE_df$log2fc
            return(out_df)
        }
    )

    names(object@data[["pairwise_DAR"]]) <- pair_names
    return(object)
})
