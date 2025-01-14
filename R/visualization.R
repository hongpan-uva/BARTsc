#' Dot plot visualization of cross-cell type analysis
#'
#' @param object A bartsc object.
#' @param mod One of "RNA", "ATAC" and "bimodal".
#' @param tf The transcription factor or chromatin regulator to inspect.
#' @param cell_types_used Cell types to use, must be a subset of the object@meta$cell_types_used.
#' @param max_dot_size Maximum size of dots.
#' @param fontsize Size of font.
#'
#' @return A plot object.
#'
#' @import ggplot2
#'
#' @export
dot_plot <- function(object, mod, tf, cell_types_used = NULL, max_dot_size = 20, fontsize = 18) {
    if (is.null(cell_types_used)) {
        cell_types_used <- object@meta$cell_types_used
    }

    if (!all(cell_types_used %in% object@meta$cell_types_used)) {
        stop("Input cell types don't match the used cell types in meta data")
    }

    cds_df <- object@resultsCrossCellType[[mod]]$CDS

    if (is.null(cds_df)) {
        stop("Queried data does not exist")
    }

    tf_cds_df <- cds_df[which(cds_df$TF == tf), ]

    if (nrow(tf_cds_df) == 0) {
        stop("Queried tf is not in the database")
    }

    plot_df <- expand.grid(rows = cell_types_used, cols = cell_types_used)
    plot_df$nlog10pvalue <- NA
    plot_df$dev <- NA
    for (i in 1:nrow(plot_df)) {
        rname <- as.character(plot_df$rows[i])
        cname <- as.character(plot_df$cols[i])
        pmtx <- object@resultsCrossCellType[[mod]]$wilcox_sign_test_pvalue[[tf]]
        dmtx <- object@resultsCrossCellType[[mod]]$deviation[[tf]]
        plot_df$nlog10pvalue[i] <- -log10(pmtx[rname, cname])
        plot_df$dev[i] <- dmtx[rname, cname]
    }

    plot_df$nlog10pvalue[which(plot_df$nlog10pvalue > 4)] <- 4

    # dot radius by -log10(pvalue)
    p <- ggplot(plot_df, aes(x = cols, y = rows, colour = dev, size = nlog10pvalue)) +
        scale_colour_gradient2(low = "#01665e", high = "#8c510a", mid = 0, name = "Deviation Ratio") +
        scale_size_area(max_size = max_dot_size, limits = c(0, 3), breaks = c(-log10(0.1), -log10(0.05), -log10(0.01), -log10(0.001)), oob = scales::squish, labels = c("0.1", "0.05", "0.01", "<0.001"), name = "P-value") +
        geom_point() +
        ggtitle(paste0(tf, " (Nprofile = ", tf_cds_df$n_profiles, ")")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        scale_y_discrete(guide = guide_axis(angle = 45)) +
        theme_bw() +
        theme(
            text = element_text(size = fontsize),
            axis.title = element_blank(),
            aspect.ratio = 1
        )

    return(p)
}


#' Heatplot visualization of cross-cell type analysis, colored by deviation ratio
#'
#' @param object A bartsc object.
#' @param mod One of "RNA", "ATAC" and "bimodal".
#' @param tf The transcription factor or chromatin regulator to inspect.
#' @param cell_types_used Cell types to use, must be a subset of the object@meta$cell_types_used.
#' @param value.display Whether to display the deviation ratio in heatmap tiles.
#' @param decimal The number of decimal places to retain for the displayed values. Default it 2.
#' @param fontsize size of font
#' @param tile_fontsize font size of displayed values in heatmap tiles.
#'
#' @return A plot object.
#'
#' @import ggplot2
#'
#' @export
deviation_heatmap <- function(
    object, mod, tf, cell_types_used = NULL,
    value.display = TRUE, decimal = 2, fontsize = 18, tile_fontsize = 7) {
    if (is.null(cell_types_used)) {
        cell_types_used <- object@meta$cell_types_used
    }

    if (!all(cell_types_used %in% object@meta$cell_types_used)) {
        stop("Input cell types don't match the used cell types in meta data")
    }

    cds_df <- object@resultsCrossCellType[[mod]]$CDS

    if (is.null(cds_df)) {
        stop("Queried data does not exist")
    }

    tf_cds_df <- cds_df[which(cds_df$TF == tf), ]

    if (nrow(tf_cds_df) == 0) {
        stop("Queried tf is not in the database")
    }

    plot_mtx <- object@resultsCrossCellType[[mod]]$deviation[[tf]][cell_types_used, cell_types_used]
    plot_df <- expand.grid(rows = rownames(plot_mtx), cols = colnames(plot_mtx))
    plot_df$dev_ratio <- 0
    for (i in 1:nrow(plot_df)) {
        row_ct <- as.character(plot_df[i, "rows"])
        col_ct <- as.character(plot_df[i, "cols"])
        plot_df$dev_ratio[i] <- plot_mtx[row_ct, col_ct]
    }
    plot_df$keep <- round(plot_df$dev_ratio, decimal)

    # heatmap
    p <- ggplot(plot_df, aes(x = cols, y = rows, fill = dev_ratio)) +
        scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = 0, name = "Deviation Ratio") +
        geom_tile(width = .9, height = .9) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        scale_y_discrete(guide = guide_axis(angle = 45)) +
        ggtitle(paste0(tf, " (Nprofle = ", tf_cds_df$n_profiles, ")")) +
        theme_bw() +
        theme(
            # panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = fontsize),
            axis.title = element_blank(),
            aspect.ratio = 1
        )

    if (value.display == TRUE) {
        p <- p + geom_text(aes(x = cols, y = rows, label = keep), color = "black", size = tile_fontsize)
    }

    return(p)
}

#' Make a 3d scatter plot for key regulators identification
#'
#' BARTsc identifies key regulators for a given cell types by comprehensively
#' considering three factors: 1) Contribution to signature expression, 2) rela-
#' tive activity across all the cell types and 3) the uniqueness of a given
#' transcription regulator. These three factors are quantified by signature
#' score, MDR (mean deviation ratio) and dMDR (differential MDR) respectively.
#' key_regulator_scatter() make a 3d scatter plot accordingly and color key
#' regulators passed p value threshold.
#'
#'
#' @param object A bartsc object.
#' @param mod One of "RNA", "ATAC" and "bimodal".
#' @param cell_type Visualize the selected cell type, must be one of the cell
#' type in object@meta$cell_types_used.
#' @param tfs_labeled Transcription regulators to be labelled on the plot.
#' @param pval.thr Threshold of the final p value.
#'
#' @return A plot.
#'
#' @import viridis
#' @import scatterplot3d
#'
#' @export
#'
key_regulator_scatter <- function(
    object, mod, cell_type, tfs_labeled = NULL, pval.thr = 0.05) {
    plot_df <- get_result(object, analysis = "Key regs ident", mod = mod)[[cell_type]]

    plot_df$role <- "other"
    if (!is.null(tfs_labeled)) {
        plot_df$role[which(plot_df$TF %in% tfs_labeled)] <- "labeled"
        plot_df$role <- factor(plot_df$role, levels = c("labeled", "other"))
    }

    n_keep <- length(which(plot_df$final_pvalue < pval.thr))
    print(paste0(n_keep, " significant key regulators were identified"))

    plot_df$label <- paste0(plot_df$final_rank, ":", plot_df$TF)
    plot_df$final_rank[which(plot_df$final_rank > n_keep)] <- NA # use top N as labeled TF

    plot_df$role2 <- rep("other", nrow(plot_df))
    plot_df$role2[which(plot_df$pval.thr < 0.05)] <- "key regulators"
    plot_df$role2 <- factor(plot_df$role2, levels = c("key regulators", "other"))

    label_df <- plot_df[which(plot_df$role == "labeled"), ]

    pallete <- rev(viridis::plasma(50)) # color pallete to use
    rank_vec <- plot_df$final_rank
    color_vec <- plot_df$final_rank
    color_boundaries <- c(min(na.omit(plot_df$final_rank)), max(na.omit(plot_df$final_rank)))
    breaks <- seq(color_boundaries[1], color_boundaries[2], length.out = 51) # color breaks

    for (i in 1:length(pallete)) {
        lower <- breaks[i]
        upper <- breaks[i + 1]
        for (j in 1:length(rank_vec)) {
            if (is.na(rank_vec[j])) {
                next()
            }
            if (rank_vec[j] >= lower && rank_vec[j] < upper) {
                color_vec[j] <- pallete[i]
            }
        }
    }
    color_vec[which(color_vec == color_boundaries[2])] <- pallete[50]
    color_vec[which(is.na(color_vec))] <- "grey"

    p <- scatterplot3d::scatterplot3d(
        x = plot_df$MDR, y = plot_df$signature_score, z = plot_df$dMDR,
        xlab = "MDR", ylab = "signature_score", zlab = "dMDR",
        pch = 16, color = color_vec,
        cex.symbols = 1.5, cex.axis = 1.5, cex.lab = 1.5,
        main = cell_type
    )

    label.coords <- p$xyz.convert(label_df$MDR, label_df$signature_score, label_df$dMDR)

    text(label.coords$x,
        label.coords$y,
        labels = label_df$label,
        cex = 1.5,
        pos = 4
    )
}
