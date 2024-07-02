setGeneric("set_used_cell_types", function(
    object,
    values) {
    standardGeneric("set_used_cell_types")
})

setMethod("set_used_cell_types", "bartsc", function(
    object,
    values) {
    if (!all(values %in% unique(object@meta$label))) {
        warning("Input cell types to use don't match exising cell types labels\n  Check exising cell types labels with object@meta$label")
    }
    object@meta$cell_types_used <- values
    return(object)
})

setGeneric("get_used_cell_types", function(
    object) {
    standardGeneric("get_used_cell_types")
})

setMethod("get_used_cell_types", "bartsc", function(
    object) {
    return(object@meta$cell_types_used)
})

setGeneric("set_active_mod", function(
    object,
    mod) {
    standardGeneric("set_active_mod")
})

setMethod("set_active_mod", "bartsc", function(
    object,
    mod) {
    if (missing(mod)) {
        stop("parameter mod can't be empty.")
    }
    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("active modality should be one of RNA, ATAC and bimodal")
    }

    object@meta$active_mod <- mod

    return(object)
})

setGeneric("get_active_mod", function(
    object) {
    standardGeneric("get_active_mod")
})

setMethod("get_active_mod", "bartsc", function(
    object) {
    return(object@meta$active_mod)
})

#' Get the analysis results for cell type signature analysis or cross-cell-
#' type analysis
#'
#' @param object A bartsc object.
#' @param analysis name of analysis to check, one of "cell type signature",
#' "cross-cell-type", "Key regs ident"
#' @param mod one of "RNA", "ATAC" and "bimodal"
#'
#' @export
#'
setGeneric("get_result", function(
    object,
    analysis = "cell type signature",
    mod) {
    standardGeneric("get_result")
})

setMethod("get_result", "bartsc", function(
    object,
    analysis = "cell type signature",
    mod) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (analysis == "cell type signature") {
        if (mod == "bimodal") {
            result <- list()
            for (ct in object@meta$cell_types_used) {
                result[[ct]] <- get_bart_result(object@resultsSignature$bimodal[[ct]], "bimodal")
            }
        } else if (mod == "RNA") {
            result <- list()
            for (ct in object@meta$cell_types_used) {
                result[[ct]] <- get_bart_result(object@resultsSignature$RNA[[ct]], "geneset")
            }
        } else if (mod == "ATAC") {
            result <- list()
            for (ct in object@meta$cell_types_used) {
                result[[ct]] <- get_bart_result(object@resultsSignature$ATAC[[ct]], "region")
            }
        }
    } else if (analysis == "cross-cell-type") {
        if (length(object@resultsCrossCellType) == 0) {
            stop("Queried data does not exist")
        }
        result <- object@resultsCrossCellType[[mod]]$deviation
    } else if (analysis == "Key regs ident") {
        if (length(object@resultsKeyRegsIdent) == 0) {
            stop("Queried data does not exist")
        }
        result <- object@resultsKeyRegsIdent[[mod]]
    }
    return(result)
})

setGeneric("get_TF_cluster", function(
    object,
    mod,
    TF) {
    standardGeneric("get_TF_cluster")
})

setMethod("get_TF_cluster", "bartsc", function(
    object,
    mod,
    TF) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]])) {
        stop("There is no existing hierarchical clustering model, run h_clustering() first.")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters"]])) {
        stop("There is no existing cut tree, run cut_tree() first.")
    }

    return(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters"]][[TF]])
})

setGeneric("get_cluster_df", function(
    object,
    mod,
    cluster) {
    standardGeneric("get_cluster_df")
})

setMethod("get_cluster_df", "bartsc", function(
    object,
    mod,
    cluster) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]])) {
        stop("There is no existing hierarchical clustering model, run h_clustering() first.")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters"]])) {
        stop("There is no existing cut tree, run cut_tree() first.")
    }

    if (missing(cluster)) {
        return(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters_df"]])
    } else {
        return(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters_df"]][[cluster]])
    }
})
