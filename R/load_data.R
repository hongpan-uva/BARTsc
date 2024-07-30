#' @export
get_RNA_cnt_from_Seurat <- function(seurat_project) {
    return(seurat_project@assays$RNA@counts)
}

#' @export
get_ATAC_cnt_from_ArchR <- function(ArchR_project, convert_name = FALSE) {
    if (!PackageCheck("ArchR")) {
        stop("package ArchR uninstalled")
    }

    peak_rse <- getMatrixFromProject(ArchRProj = ArchR_project, useMatrix = "PeakMatrix")
    peak_mtx <- assays(peak_rse)[[1]]
    rownames(peak_mtx) <- 1:nrow(peak_mtx) # give row names

    if (convert_name == TRUE) {
        cnames <- gsub("-", ".", colnames(peak_mtx)) # rename cells
        cnames <- sapply(cnames, function(x) {
            return(strsplit(x, "#")[[1]][2])
        }) # rename cells
        colnames(peak_mtx) <- cnames
    }

    return(peak_mtx)
}

#' @export
get_peaks_from_ArchR <- function(ArchR_project) {
    if (!PackageCheck("ArchR")) {
        stop("package ArchR uninstalled")
    }

    peaks_gr <- getPeakSet(ArchRProj = ArchR_project)
    peaks <- data.frame(peaks_gr)
    peaks <- cbind(peaks[, 1:3], name = ".", score = ".", strand = ".")
    return(peaks)
}

#' @export
get_label_from_Seurat <- function(seurat_project, col) {
    label <- proj@meta.data[[col]]
    names(label) <- rownames(proj@meta.data)
    return(label)
}

#' @export
get_label_from_ArchR <- function(ArchR_project, col, convert_name = FALSE) {
    Cdata <- getCellColData(ArchR_project)
    label <- Cdata[, col]

    if (convert_name == TRUE) {
        cnames <- gsub("-", ".", rownames(Cdata)) # rename cells
        cnames <- sapply(cnames, function(x) {
            return(strsplit(x, "#")[[1]][2])
        }) # rename cells
        names(label) <- cnames
    } else {
        names(label) <- rownames(Cdata)
    }

    return(label)
}
