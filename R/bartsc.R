#' An S4 class to store raw data, intermediate data and outcomes of BARTsc analysis
#'
#' @slot meta, param, data, Args, intermediate, result
setClass("bartsc",
    slots = c(
        meta = "list",
        param = "list",
        assays = "list",
        data = "list",
        resultsSignature = "list",
        resultsCrossCellType = "list",
        resultsKeyRegsIdent = "list"
    ),
    prototype = list(
        meta = list(),
        param = list(),
        assays = list(),
        data = list(),
        resultsSignature = list(),
        resultsCrossCellType = list(),
        resultsKeyRegsIdent = list()
    )
)

#' Create a bartsc object.
#'
#' bartsc object is used to store input data, processed data and results of
#' BARTsc analyses. With required input data, this function creates a bartsc
#' object of scRNA-seq, scATAC-seq or the bimodal sequencing of transcriptiome
#' and chromatin accessibility (single-cell multiomics).
#'
#' @param name Name to use for the object.
#' @param genome Genome to use, hg38 or mm10.
#' @param label A vector of cell types, with cell ids as names.
#' @param cell_types_used A vector of cell types to use. cell types must exist
#' in `label`.
#' @param RNA_cnt_matrix A gene by cell read count matrix
#' @param RNA_norm_matrix A gene by cell normalized gene expression matrix.
#' Must be a `dgCMatrix` object.
#' @param ATAC_cnt_matrix A peak by cell count matrix. Must be a `dgCMatrix`
#' object.
#' @param ATAC_norm_matrix  A peak by cell normalized chromatin
#' accessibility matrix. Must be a `dgCMatrix` object.
#' @param peaks A dataframe of input peaks.
#' @param gene_mode_param A list of parameters for gene mode.
#' @param region_mode_param A list of parameters list for region mode.
#' @param bimodal_mode_param A list of parameters list for bimodal mode.
#'
#' @return a bartsc object
#'
#' @export
#'
bartsc <- function(name, genome, label, cell_types_used = NULL,
                   RNA_cnt_matrix = NULL, RNA_norm_matrix = NULL,
                   ATAC_cnt_matrix = NULL, ATAC_norm_matrix = NULL,
                   peaks = NULL,
                   gene_mode_param = list(binsize = 1000),
                   region_mode_param = list(
                       binsize = 50, scorecol = 5
                   ),
                   bimodal_mode_param = list(binsize = 50)) {
    # input validation
    stopifnot(exprs = {
        !missing(name)
        !missing(genome)
        genome %in% c("hg38", "mm10")
        !missing(label)
    })

    if ((!is.null(ATAC_cnt_matrix) || !is.null(ATAC_norm_matrix)) != !is.null(peaks)) {
        stop("count matrix or normalized matrix of ATAC-seq must be input with peaks")
    }

    active_mod <- NULL
    rna_input <- FALSE
    atac_input <- FALSE

    if (!is.null(RNA_cnt_matrix) || !is.null(RNA_norm_matrix)) {
        rna_input <- TRUE
    }
    if (
        (!is.null(ATAC_cnt_matrix) || !is.null(ATAC_norm_matrix)) &&
            !is.null(peaks)
    ) {
        atac_input <- TRUE
    }

    # check if cell names match
    if (rna_input == TRUE && atac_input == TRUE) {
        if (!all(colnames(RNA_cnt_matrix) %in% colnames(ATAC_cnt_matrix))) {
            stop("cell names of RNA modality and ATAC modality do not match")
        }
    }

    object <- new("bartsc")
    object@meta$name <- name
    object@meta$genome <- genome

    if (is.null(cell_types_used)) {
        cell_types_used <- unique(as.character(label))
    }
    object@meta$cell_types_used <- cell_types_used

    if (class(label) != "factor") {
        object@meta$label <- factor(label, levels = cell_types_used)
    } else {
        object@meta$label <- label
    }

    # set geneset mode data and paramter
    if (!is.null(RNA_cnt_matrix)) {
        object@assays[["RNA"]][["counts"]] <- RNA_cnt_matrix
    } else if (!is.null(RNA_norm_matrix)) {
        object@assays[["RNA"]][["normalized"]] <- RNA_norm_matrix
    }

    # set region mode data
    if (!is.null(ATAC_cnt_matrix) && !is.null(peaks)) {
        object@assays[["ATAC"]] <- list()
        object@assays[["ATAC"]][["peaks"]] <- peaks
        rownames(ATAC_cnt_matrix) <- 1:nrow(ATAC_cnt_matrix) # coerce rownames
        object@assays[["ATAC"]][["counts"]] <- ATAC_cnt_matrix
    } else if (!is.null(ATAC_norm_matrix) && !is.null(peaks)) {
        object@assays[["ATAC"]] <- list()
        object@assays[["ATAC"]][["peaks"]] <- peaks
        object@assays[["ATAC"]][["normalized"]] <- ATAC_norm_matrix
    }

    # set paramters
    object@param[["gene_mode_param"]] <- gene_mode_param
    object@param[["region_mode_param"]] <- region_mode_param
    object@param[["bimodal_mode_param"]] <- bimodal_mode_param

    # set default active modality
    if (rna_input == TRUE && atac_input == TRUE) {
        active_mod <- "bimodal"
    } else if (rna_input == TRUE) {
        active_mod <- "RNA"
    } else if (atac_input == TRUE) {
        active_mod <- "ATAC"
    }

    object@meta$active_mod <- active_mod

    return(object)
}

#' Normalize the scRNA-seq count matrix in a bartsc object.
#'
#' Feature counts for each cell are divided by the total
#' counts for that cell and multiplied by 10,000. This
#' is then natural-log transformed using ‘log1p’.
#'
#' @param object A bartsc object.
#'
#' @return A bartsc object.
#'
#' @importFrom methods as
#'
#' @export
#'
setGeneric("normalize_RNA", function(
    object) {
    standardGeneric("normalize_RNA")
})

setMethod("normalize_RNA", "bartsc", function(
    object) {
    if (is.null(object@assays[["RNA"]][["counts"]])) {
        stop("no RNA counts matrix in bartsc object")
    }
    scale.factor <- 10000
    # perform log normalization
    mtx_by_cnt <- apply(object@assays[["RNA"]][["counts"]], 2, function(x) {
        read_cnt <- sum(x)
        return(x / read_cnt)
    })
    mtx_norm <- log1p(scale.factor * mtx_by_cnt)
    mtx_norm <- as(mtx_norm, "dgCMatrix")
    object@assays[["RNA"]][["normalized"]] <- mtx_norm
    return(object)
})

#' Normalize the scATAC-seq count matrix in a bartsc object.
#'
#' Use TF-IDF to normalize scATAC-seq data. scale factor is set as 1e4.
#'
#' @param object A bartsc object.
#' @param method Which method to use for calculation. To check details,
#' use help(RunTFIDF).
#' @param scale.factor Which scale factor to use. Default is 10000.
#'
#' @return A bartsc object.
#'
#' @importFrom methods as
#'
#' @export
#'
setGeneric("normalize_ATAC", function(
    object, method.use = 1, scale.factor = 1e4) {
    standardGeneric("normalize_ATAC")
})

setMethod("normalize_ATAC", "bartsc", function(
    object, method.use = 1, scale.factor = 1e4) {
    if (is.null(object@assays[["ATAC"]][["counts"]])) {
        stop("no ATAC counts matrix in bartsc object")
    }
    object@assays[["ATAC"]][["normalized"]] <- RunTFIDF(object@assays[["ATAC"]][["counts"]], method = method.use, scale.factor = scale.factor, idf = NULL, verbose = TRUE)
    return(object)
})

#' Run cell type signature analysis on RNA assay
#'
#' @param object A bartsc object.
#'
#' @export
#'
setGeneric("run_signature_RNA", function(object) standardGeneric("run_signature_RNA"))

setMethod("run_signature_RNA", "bartsc", function(object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))
    object@resultsSignature[["RNA"]] <- list()
    bart_list <- mclapply(
        X = object@meta$cell_types_used,
        FUN = function(x) {
            bart_obj <- bart(x, object@meta$genome,
                gene_data = object@data$signature_genes[[x]],
                gene_mode_param = object@param$gene_mode_param
            )
            bart_obj <- runBartGeneSet(bart_obj)
            return(bart_obj)
        }
    )
    names(bart_list) <- object@meta$cell_types_used
    object@resultsSignature[["RNA"]] <- bart_list
    return(object)
})

#' Run cell type signature analysis on ATAC assay
#'
#' @param object A bartsc object.
#'
#' @export
#'
setGeneric("run_signature_ATAC", function(object) standardGeneric("run_signature_ATAC"))

setMethod("run_signature_ATAC", "bartsc", function(object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))
    object@resultsSignature[["ATAC"]] <- list()
    bart_list <- mclapply(
        X = object@meta$cell_types_used,
        FUN = function(x) {
            bart_obj <- bart(x, object@meta$genome,
                region_data = object@data$signature_peaks[[x]],
                region_mode_param = object@param$region_mode_param
            )
            bart_obj <- runBartRegion(bart_obj)
            return(bart_obj)
        }
    )
    names(bart_list) <- object@meta$cell_types_used
    object@resultsSignature[["ATAC"]] <- bart_list
    return(object)
})

#' Run cell type signature analysis on RNA assay and ATAC assay with bimodal integration method
#'
#' @param object A bartsc object.
#'
#' @export
#'
setGeneric("run_signature_bimodal", function(object) standardGeneric("run_signature_bimodal"))

setMethod("run_signature_bimodal", "bartsc", function(object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))
    object@resultsSignature[["bimodal"]] <- list()
    bart_list <- mclapply(
        X = object@meta$cell_types_used,
        FUN = function(x) {
            bart_obj <- bart(x, object@meta$genome,
                gene_data = object@data$signature_genes[[x]],
                region_data = object@data$signature_peaks[[x]],
                bimodal_mode_param = object@param$bimodal_mode_param
            )
            bart_obj <- runBartBimodal(bart_obj)
            return(bart_obj)
        }
    )
    names(bart_list) <- object@meta$cell_types_used
    object@resultsSignature[["bimodal"]] <- bart_list
    return(object)
})

#' calculate aucs for cross-cell-type analysis of RNA
#'
#' @param object a bart object
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
#'
#' @export
#'
setGeneric("calc_crossCT_auc_RNA", function(
    object) {
    standardGeneric("calc_crossCT_auc_RNA")
})

setMethod("calc_crossCT_auc_RNA", "bartsc", function(
    object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))

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

    object@data[["RNA_pairwise_AUC"]] <- mclapply(
        X = pairs,
        FUN = function(x) {
            message(x)
            bart_obj <- bart(paste0(x[1], "::", x[2]),
                genome = object@meta$genome,
                gene_data = object@data[["pairwise_DEG"]][[paste0(x[1], "::", x[2])]],
                gene_mode_param = object@param$gene_mode_param
            )
            bart_obj <- runBartGeneSet(bart_obj)
            auc_df <- bart_obj@result[["geneset"]][["auc"]]
            return(auc_df)
        }
    )

    names(object@data[["RNA_pairwise_AUC"]]) <- pair_names
    return(object)
})

#' calculate aucs for cross-cell-type analysis of bimodal
#'
#' @param object a bart object
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
#'
#' @export
#'
setGeneric("calc_crossCT_auc_bimodal", function(
    object) {
    standardGeneric("calc_crossCT_auc_bimodal")
})

setMethod("calc_crossCT_auc_bimodal", "bartsc", function(
    object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))

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

    object@data[["bimodal_pairwise_AUC"]] <- mclapply(
        X = pairs,
        FUN = function(x) {
            message(x)
            bart_obj <- bart(paste0(x[1], "::", x[2]),
                genome = object@meta$genome,
                gene_data = object@data[["pairwise_DEG"]][[paste0(x[1], "::", x[2])]],
                region_data = object@data[["pairwise_DAR"]][[paste0(x[1], "::", x[2])]],
                bimodal_mode_param = object@param$bimodal_mode_param
            )
            bart_obj <- runBartBimodal(bart_obj)
            auc_df <- bart_obj@result[["bimodal"]][["auc"]]
            return(auc_df)
        }
    )

    names(object@data[["bimodal_pairwise_AUC"]]) <- pair_names
    return(object)
})

#' calculate aucs for cross-cell-type analysis of ATAC
#'
#' @param object a bart object
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
#'
#' @export
#'
setGeneric("calc_crossCT_auc_ATAC", function(
    object) {
    standardGeneric("calc_crossCT_auc_ATAC")
})

setMethod("calc_crossCT_auc_ATAC", "bartsc", function(
    object) {
    message(paste(getOption("mc.cores", 2L), "cores to use"))

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

    object@data[["ATAC_pairwise_AUC"]] <- mclapply(
        X = pairs,
        FUN = function(x) {
            message(x)
            bart_obj <- bart(paste0(x[1], "::", x[2]),
                genome = object@meta$genome,
                region_data = object@data[["pairwise_DAR"]][[paste0(x[1], "::", x[2])]],
                region_mode_param = object@param$region_mode_param
            )
            bart_obj <- runBartRegion(bart_obj)
            auc_df <- bart_obj@result[["region"]][["auc"]]
            return(auc_df)
        }
    )

    names(object@data[["ATAC_pairwise_AUC"]]) <- pair_names
    return(object)
})

#' create bart objects and run BART for cross-cell-type test
#'
#' @param object a bart object
#' @param mod one of "RNA", "ATAC" and "bimodal"
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
#'
#' @export
#'
setGeneric("crossCT_test", function(
    object,
    mod) {
    standardGeneric("crossCT_test")
})

setMethod("crossCT_test", "bartsc", function(
    object,
    mod = NULL) {
    if (is.null(mod) || !(mod %in% c("RNA", "ATAC", "bimodal"))) {
        stop('mod must be one of "RNA", "ATAC" and "bimodal"')
    }

    if (mod == "RNA") {
        auc_title <- "RNA_pairwise_AUC"
    } else if (mod == "ATAC") {
        auc_title <- "ATAC_pairwise_AUC"
    } else if (mod == "bimodal") {
        auc_title <- "bimodal_pairwise_AUC"
    }

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

    # prepare auc data frames
    auc_df <- as.data.frame(matrix(nrow = 0, ncol = 4))

    for (i in pair_names) {
        auc <- object@data[[auc_title]][[i]]
        auc$rank <- rank(auc$AUC) # get rank, greater number for greater auc
        auc$TF <- sapply(auc$ChIP_seq, function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        auc$cluster_pair <- i
        auc$zscore <- (auc$AUC - mean(auc$AUC)) / sd(auc$AUC)
        auc_df <- rbind(auc_df, auc)
    }

    tf_list <- unique(auc_df$TF)

    # pairwise wilcoxon test
    wcx_stat_list <- list()
    wcx_pvalue_list <- list()
    ratio_mtx_list <- list()
    dev_mtx_list <- list()
    avg_abs_dev_list <- list()
    cds_df <- data.frame(matrix(nrow = 0, ncol = 4)) # Average Differential Score

    for (tf in tf_list) {
        tmp_df <- auc_df[which(auc_df$TF == tf), ]

        wcx_stat <- matrix(nrow = length(cell_types_used), ncol = length(cell_types_used))
        rownames(wcx_stat) <- cell_types_used
        colnames(wcx_stat) <- cell_types_used

        wcx_pvalue <- wcx_stat

        n_profiles <- length(unique(tmp_df$ChIP_seq))
        max_rank_sum <- n_profiles * (n_profiles + 1) / 2

        # pairwise.wilcox.test(tmp_df$AUC, tmp_df$cluster, alternative = "less", paired = TRUE)
        for (i in cell_types_used) {
            for (j in cell_types_used) {
                if (i == j) {
                    wcx_stat[i, j] <- max_rank_sum / 2
                    wcx_pvalue[i, j] <- 1
                    next
                }
                # vec1 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                # vec1 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                vec1 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                vec2 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                test.res <- wilcox.test(vec1, vec2,
                    paired = TRUE, alternative = "two.sided"
                )
                wcx_stat[i, j] <- test.res$statistic
                wcx_pvalue[i, j] <- test.res$p.value
            }
        }

        ratio_mtx <- wcx_stat / max_rank_sum
        deviation_mtx <- ratio_mtx * 2 - 1
        abs_deviation_mtx <- abs(deviation_mtx)
        avg_abs_dev <- mean(abs_deviation_mtx[upper.tri(abs_deviation_mtx)])
        n_pairs <- length(wcx_pvalue[which(wcx_pvalue < 0.05)])

        wcx_stat_list[[tf]] <- wcx_stat
        wcx_pvalue_list[[tf]] <- wcx_pvalue
        ratio_mtx_list[[tf]] <- ratio_mtx
        dev_mtx_list[[tf]] <- deviation_mtx
        avg_abs_dev_list[[tf]] <- avg_abs_dev
        cds_df <- rbind(cds_df, c(tf, avg_abs_dev, n_pairs, n_profiles)) # cds for Comprehensive Deviation Score
    }

    colnames(cds_df) <- c("TF", "CDS", "n_sgfnt_pairs", "n_profiles")
    rownames(cds_df) <- 1:nrow(cds_df)
    cds_df$CDS <- as.numeric(cds_df$CDS)
    cds_df$n_sgfnt_pairs <- as.integer(cds_df$n_sgfnt_pairs)
    cds_df$n_profiles <- as.integer(cds_df$n_profiles)
    cds_df <- cds_df[order(cds_df$n_sgfnt_pairs, cds_df$CDS, decreasing = TRUE), ]

    object@resultsCrossCellType[[mod]] <- list(
        CDS = cds_df, deviation = dev_mtx_list, wilcox_sign_test_pvalue = wcx_pvalue_list
    )

    return(object)
})

#' create bart objects and run BART for cross-cell-type test
#'
#' @param object a bart object
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
setGeneric("crossCT_test_RNA", function(
    object) {
    standardGeneric("crossCT_test_RNA")
})

setMethod("crossCT_test_RNA", "bartsc", function(
    object) {
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

    # prepare auc data frames
    auc_df <- as.data.frame(matrix(nrow = 0, ncol = 4))

    for (i in pair_names) {
        auc <- object@data$RNA_pairwise_AUC[[i]]
        auc$rank <- rank(auc$AUC) # get rank, greater number for greater auc
        auc$TF <- sapply(auc$ChIP_seq, function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        auc$cluster_pair <- i
        auc$zscore <- (auc$AUC - mean(auc$AUC)) / sd(auc$AUC)
        auc_df <- rbind(auc_df, auc)
    }

    tf_list <- unique(auc_df$TF)

    # pairwise wilcoxon test
    wcx_stat_list <- list()
    wcx_pvalue_list <- list()
    ratio_mtx_list <- list()
    dev_mtx_list <- list()
    avg_abs_dev_list <- list()
    cds_df <- data.frame(matrix(nrow = 0, ncol = 4)) # Average Differential Score

    for (tf in tf_list) {
        tmp_df <- auc_df[which(auc_df$TF == tf), ]

        wcx_stat <- matrix(nrow = length(cell_types_used), ncol = length(cell_types_used))
        rownames(wcx_stat) <- cell_types_used
        colnames(wcx_stat) <- cell_types_used

        wcx_pvalue <- wcx_stat

        n_profiles <- length(unique(tmp_df$ChIP_seq))
        max_rank_sum <- n_profiles * (n_profiles + 1) / 2

        # pairwise.wilcox.test(tmp_df$AUC, tmp_df$cluster, alternative = "less", paired = TRUE)
        for (i in cell_types_used) {
            for (j in cell_types_used) {
                if (i == j) {
                    wcx_stat[i, j] <- max_rank_sum / 2
                    wcx_pvalue[i, j] <- 1
                    next
                }
                # vec1 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                # vec1 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                vec1 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                vec2 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                test.res <- wilcox.test(vec1, vec2,
                    paired = TRUE, alternative = "two.sided"
                )
                wcx_stat[i, j] <- test.res$statistic
                wcx_pvalue[i, j] <- test.res$p.value
            }
        }

        ratio_mtx <- wcx_stat / max_rank_sum
        deviation_mtx <- ratio_mtx * 2 - 1
        abs_deviation_mtx <- abs(deviation_mtx)
        avg_abs_dev <- mean(abs_deviation_mtx[upper.tri(abs_deviation_mtx)])
        n_pairs <- length(wcx_pvalue[which(wcx_pvalue < 0.05)])

        wcx_stat_list[[tf]] <- wcx_stat
        wcx_pvalue_list[[tf]] <- wcx_pvalue
        ratio_mtx_list[[tf]] <- ratio_mtx
        dev_mtx_list[[tf]] <- deviation_mtx
        avg_abs_dev_list[[tf]] <- avg_abs_dev
        cds_df <- rbind(cds_df, c(tf, avg_abs_dev, n_pairs, n_profiles)) # cds for Comprehensive Deviation Score
    }

    colnames(cds_df) <- c("TF", "CDS", "n_sgfnt_pairs", "n_profiles")
    rownames(cds_df) <- 1:nrow(cds_df)
    cds_df$CDS <- as.numeric(cds_df$CDS)
    cds_df$n_sgfnt_pairs <- as.integer(cds_df$n_sgfnt_pairs)
    cds_df$n_profiles <- as.integer(cds_df$n_profiles)
    cds_df <- cds_df[order(cds_df$n_sgfnt_pairs, cds_df$CDS, decreasing = TRUE), ]

    object@resultsCrossCellType[["RNA"]] <- list(
        CDS = cds_df, deviation = dev_mtx_list, wilcox_sign_test_pvalue = wcx_pvalue_list
    )

    return(object)
})

#' create bart objects and run BART for cross-cell-type test
#'
#' @param object a bart object
#'
#' @importFrom combinat combn
#'
#' @return a bartsc object
setGeneric("crossCT_test_bimodal", function(
    object) {
    standardGeneric("crossCT_test_bimodal")
})

setMethod("crossCT_test_bimodal", "bartsc", function(
    object) {
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

    # prepare auc data frames
    auc_df <- as.data.frame(matrix(nrow = 0, ncol = 4))

    for (i in pair_names) {
        auc <- object@data$bimodal_pairwise_AUC[[i]]
        auc$rank <- rank(auc$AUC) # get rank, greater number for greater auc
        auc$TF <- sapply(auc$ChIP_seq, function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        auc$cluster_pair <- i
        auc$zscore <- (auc$AUC - mean(auc$AUC)) / sd(auc$AUC)
        auc_df <- rbind(auc_df, auc)
    }

    tf_list <- unique(auc_df$TF)

    # pairwise wilcoxon test
    wcx_stat_list <- list()
    wcx_pvalue_list <- list()
    ratio_mtx_list <- list()
    dev_mtx_list <- list()
    avg_abs_dev_list <- list()
    cds_df <- data.frame(matrix(nrow = 0, ncol = 4)) # Average Differential Score

    for (tf in tf_list) {
        tmp_df <- auc_df[which(auc_df$TF == tf), ]

        wcx_stat <- matrix(nrow = length(cell_types_used), ncol = length(cell_types_used))
        rownames(wcx_stat) <- cell_types_used
        colnames(wcx_stat) <- cell_types_used

        wcx_pvalue <- wcx_stat

        n_profiles <- length(unique(tmp_df$ChIP_seq))
        max_rank_sum <- n_profiles * (n_profiles + 1) / 2

        # pairwise.wilcox.test(tmp_df$AUC, tmp_df$cluster, alternative = "less", paired = TRUE)
        for (i in cell_types_used) {
            for (j in cell_types_used) {
                if (i == j) {
                    wcx_stat[i, j] <- max_rank_sum / 2
                    wcx_pvalue[i, j] <- 1
                    next
                }
                # vec1 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                # vec1 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                vec1 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                vec2 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                test.res <- wilcox.test(vec1, vec2,
                    paired = TRUE, alternative = "two.sided"
                )
                wcx_stat[i, j] <- test.res$statistic
                wcx_pvalue[i, j] <- test.res$p.value
            }
        }

        ratio_mtx <- wcx_stat / max_rank_sum
        deviation_mtx <- ratio_mtx * 2 - 1
        abs_deviation_mtx <- abs(deviation_mtx)
        avg_abs_dev <- mean(abs_deviation_mtx[upper.tri(abs_deviation_mtx)])
        n_pairs <- length(wcx_pvalue[which(wcx_pvalue < 0.05)])

        wcx_stat_list[[tf]] <- wcx_stat
        wcx_pvalue_list[[tf]] <- wcx_pvalue
        ratio_mtx_list[[tf]] <- ratio_mtx
        dev_mtx_list[[tf]] <- deviation_mtx
        avg_abs_dev_list[[tf]] <- avg_abs_dev
        cds_df <- rbind(cds_df, c(tf, avg_abs_dev, n_pairs, n_profiles)) # cds for Comprehensive Deviation Score
    }

    colnames(cds_df) <- c("TF", "CDS", "n_sgfnt_pairs", "n_profiles")
    rownames(cds_df) <- 1:nrow(cds_df)
    cds_df$CDS <- as.numeric(cds_df$CDS)
    cds_df$n_sgfnt_pairs <- as.integer(cds_df$n_sgfnt_pairs)
    cds_df$n_profiles <- as.integer(cds_df$n_profiles)
    cds_df <- cds_df[order(cds_df$n_sgfnt_pairs, cds_df$CDS, decreasing = TRUE), ]

    object@resultsCrossCellType[["bimodal"]] <- list(
        CDS = cds_df, deviation = dev_mtx_list, wilcox_sign_test_pvalue = wcx_pvalue_list
    )

    return(object)
})

#' create bart objects and run BART for cross-cell-type test
#'
#' @param object a bart object
#'
#' @return a bartsc object
setGeneric("crossCT_test_ATAC", function(
    object) {
    standardGeneric("crossCT_test_ATAC")
})

setMethod("crossCT_test_ATAC", "bartsc", function(
    object) {
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

    # prepare auc data frames
    auc_df <- as.data.frame(matrix(nrow = 0, ncol = 4))

    for (i in pair_names) {
        auc <- object@data$ATAC_pairwise_AUC[[i]]
        auc$rank <- rank(auc$AUC) # get rank, greater number for greater auc
        auc$TF <- sapply(auc$ChIP_seq, function(x) {
            return(strsplit(x, "_")[[1]][1])
        })
        auc$cluster_pair <- i
        auc$zscore <- (auc$AUC - mean(auc$AUC)) / sd(auc$AUC)
        auc_df <- rbind(auc_df, auc)
    }

    tf_list <- unique(auc_df$TF)

    # pairwise wilcoxon test
    wcx_stat_list <- list()
    wcx_pvalue_list <- list()
    ratio_mtx_list <- list()
    dev_mtx_list <- list()
    avg_abs_dev_list <- list()
    cds_df <- data.frame(matrix(nrow = 0, ncol = 4)) # Average Differential Score

    for (tf in tf_list) {
        tmp_df <- auc_df[which(auc_df$TF == tf), ]

        wcx_stat <- matrix(nrow = length(cell_types_used), ncol = length(cell_types_used))
        rownames(wcx_stat) <- cell_types_used
        colnames(wcx_stat) <- cell_types_used

        wcx_pvalue <- wcx_stat

        n_profiles <- length(unique(tmp_df$ChIP_seq))
        max_rank_sum <- n_profiles * (n_profiles + 1) / 2

        # pairwise.wilcox.test(tmp_df$AUC, tmp_df$cluster, alternative = "less", paired = TRUE)
        for (i in cell_types_used) {
            for (j in cell_types_used) {
                if (i == j) {
                    wcx_stat[i, j] <- max_rank_sum / 2
                    wcx_pvalue[i, j] <- 1
                    next
                }
                # vec1 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$AUC[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                # vec1 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                # vec2 <- tmp_df$rank[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                vec1 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(i, "::", j))]
                vec2 <- tmp_df$zscore[which(tmp_df$cluster_pair == paste0(j, "::", i))]
                test.res <- wilcox.test(vec1, vec2,
                    paired = TRUE, alternative = "two.sided"
                )
                wcx_stat[i, j] <- test.res$statistic
                wcx_pvalue[i, j] <- test.res$p.value
            }
        }

        ratio_mtx <- wcx_stat / max_rank_sum
        deviation_mtx <- ratio_mtx * 2 - 1
        abs_deviation_mtx <- abs(deviation_mtx)
        avg_abs_dev <- mean(abs_deviation_mtx[upper.tri(abs_deviation_mtx)])
        n_pairs <- length(wcx_pvalue[which(wcx_pvalue < 0.05)])

        wcx_stat_list[[tf]] <- wcx_stat
        wcx_pvalue_list[[tf]] <- wcx_pvalue
        ratio_mtx_list[[tf]] <- ratio_mtx
        dev_mtx_list[[tf]] <- deviation_mtx
        avg_abs_dev_list[[tf]] <- avg_abs_dev
        cds_df <- rbind(cds_df, c(tf, avg_abs_dev, n_pairs, n_profiles)) # cds for Comprehensive Deviation Score
    }

    colnames(cds_df) <- c("TF", "CDS", "n_sgfnt_pairs", "n_profiles")
    rownames(cds_df) <- 1:nrow(cds_df)
    cds_df$CDS <- as.numeric(cds_df$CDS)
    cds_df$n_sgfnt_pairs <- as.integer(cds_df$n_sgfnt_pairs)
    cds_df$n_profiles <- as.integer(cds_df$n_profiles)
    cds_df <- cds_df[order(cds_df$n_sgfnt_pairs, cds_df$CDS, decreasing = TRUE), ]

    object@resultsCrossCellType[["ATAC"]] <- list(
        CDS = cds_df, deviation = dev_mtx_list, wilcox_sign_test_pvalue = wcx_pvalue_list
    )

    return(object)
})
