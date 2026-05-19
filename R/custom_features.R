#' Set custom feature lists for BARTsc analysis
#'
#' Manually assign signature genes, pairwise DEGs, signature peaks,
#' and pairwise DARs to a bartsc object. This bypasses the built-in
#' feature extraction functions.
#'
#' @param object A bartsc object.
#' @param signature_genes Named list of character vectors (gene symbols).
#'   Names must match cell types in \code{object@meta$cell_types_used}.
#' @param pairwise_DEG Named list of character vectors.
#'   Names must be ordered pairs in \code{"ct1::ct2"} format.
#' @param signature_peaks Named list of data frames (BED-like format).
#' @param pairwise_DAR Named list of data frames (BED-like format).
#' @param validate Logical. If TRUE, validates names against object metadata.
#'
#' @return A bartsc object.
#'
#' @export
set_custom_features <- function(object,
                                signature_genes = NULL,
                                pairwise_DEG = NULL,
                                signature_peaks = NULL,
                                pairwise_DAR = NULL,
                                validate = TRUE) {
    cell_types_used <- object@meta$cell_types_used

    if (!is.null(signature_genes)) {
        if (!is.list(signature_genes)) {
            stop("signature_genes must be a named list")
        }
        if (validate) {
            missing_ct <- setdiff(cell_types_used, names(signature_genes))
            extra_ct <- setdiff(names(signature_genes), cell_types_used)
            if (length(missing_ct) > 0) {
                warning(paste("Missing signature_genes for cell types:",
                              paste(missing_ct, collapse = ", ")))
            }
            if (length(extra_ct) > 0) {
                warning(paste("Extra signature_genes entries not in cell_types_used:",
                              paste(extra_ct, collapse = ", ")))
            }
        }
        object@data$signature_genes <- signature_genes
    }

    if (!is.null(pairwise_DEG)) {
        if (!is.list(pairwise_DEG)) {
            stop("pairwise_DEG must be a named list")
        }
        if (validate) {
            expected_pairs <- c()
            for (ct1 in cell_types_used) {
                for (ct2 in cell_types_used) {
                    if (ct1 != ct2) {
                        expected_pairs <- c(expected_pairs, paste0(ct1, "::", ct2))
                    }
                }
            }
            missing_pairs <- setdiff(expected_pairs, names(pairwise_DEG))
            extra_pairs <- setdiff(names(pairwise_DEG), expected_pairs)
            if (length(missing_pairs) > 0) {
                warning(paste("Missing pairwise_DEG for pairs:",
                              paste(missing_pairs, collapse = ", ")))
            }
            if (length(extra_pairs) > 0) {
                warning(paste("Extra pairwise_DEG entries:",
                              paste(extra_pairs, collapse = ", ")))
            }
        }
        object@data$pairwise_DEG <- pairwise_DEG
    }

    if (!is.null(signature_peaks)) {
        if (!is.list(signature_peaks)) {
            stop("signature_peaks must be a named list")
        }
        if (validate) {
            missing_ct <- setdiff(cell_types_used, names(signature_peaks))
            extra_ct <- setdiff(names(signature_peaks), cell_types_used)
            if (length(missing_ct) > 0) {
                warning(paste("Missing signature_peaks for cell types:",
                              paste(missing_ct, collapse = ", ")))
            }
            if (length(extra_ct) > 0) {
                warning(paste("Extra signature_peaks entries not in cell_types_used:",
                              paste(extra_ct, collapse = ", ")))
            }
        }
        object@data$signature_peaks <- signature_peaks
    }

    if (!is.null(pairwise_DAR)) {
        if (!is.list(pairwise_DAR)) {
            stop("pairwise_DAR must be a named list")
        }
        if (validate) {
            expected_pairs <- c()
            for (ct1 in cell_types_used) {
                for (ct2 in cell_types_used) {
                    if (ct1 != ct2) {
                        expected_pairs <- c(expected_pairs, paste0(ct1, "::", ct2))
                    }
                }
            }
            missing_pairs <- setdiff(expected_pairs, names(pairwise_DAR))
            extra_pairs <- setdiff(names(pairwise_DAR), expected_pairs)
            if (length(missing_pairs) > 0) {
                warning(paste("Missing pairwise_DAR for pairs:",
                              paste(missing_pairs, collapse = ", ")))
            }
            if (length(extra_pairs) > 0) {
                warning(paste("Extra pairwise_DAR entries:",
                              paste(extra_pairs, collapse = ", ")))
            }
        }
        object@data$pairwise_DAR <- pairwise_DAR
    }

    return(object)
}
