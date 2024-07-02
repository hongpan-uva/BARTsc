setGeneric("h_clustering", function(
    object,
    mod,
    min_n_profile = 3) {
    standardGeneric("h_clustering")
})

setMethod("h_clustering", "bartsc", function(
    object,
    mod,
    min_n_profile = 3) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    CDS <- object@resultsCrossCellType[[mod]]$CDS
    CDS <- CDS[which(CDS$n_profiles >= min_n_profile), ]

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

    ratio_vec_list <- list()

    dev_mtx_list <- object@resultsCrossCellType[[mod]]$deviation

    for (i in CDS$TF) {
        vec <- c()
        mtx <- dev_mtx_list[[i]]
        dims <- nrow(mtx)
        for (j in 1:dims) {
            vec <- c(vec, mtx[j, -j])
        }
        vec <- unname(vec)

        ratio_vec_list[[i]] <- vec
    }

    df <- as.data.frame(ratio_vec_list) %>% t()
    colnames(df) <- pair_names
    object@resultsCrossCellType[[mod]][["h_clustering"]][["deviation_df"]] <- df

    model <- hclust(dist(df), "average")
    object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]] <- model

    return(object)
})

setGeneric("cut_tree", function(
    object,
    mod,
    cut_height) {
    standardGeneric("cut_tree")
})

setMethod("cut_tree", "bartsc", function(
    object,
    mod,
    cut_height) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]])) {
        stop("No hierarchical clustering result. Run h_clustering(object, mod) first.")
    } else {
        model <- object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]]
    }

    if (cut_height < 0 || cut_height > max(model$height)) {
        stop(paste0("valid cut height is from 0 to ", max(model$height)))
    }

    cut_clusters <- cutree(model, h = cut_height)

    object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_height"]] <- cut_height
    object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters"]] <- cut_clusters

    CDS <- object@resultsCrossCellType[[mod]]$CDS

    cut_clusters_df <- list()
    for (i in 1:max(cut_clusters)) {
        cluster_TFs <- names(cut_clusters)[which(cut_clusters == i)] %>%
            gsub("\\.", "-", .)
        CDS_cluster <- CDS[which(CDS$TF %in% cluster_TFs), ]
        rownames(CDS_cluster) <- 1:nrow(CDS_cluster)
        cut_clusters_df[[i]] <- CDS_cluster
    }

    object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_clusters_df"]] <- cut_clusters_df

    return(object)
})

setGeneric("plot_clusters", function(
    object,
    mod,
    font_size = 11,
    cut = TRUE,
    cut_height,
    gap_size = 5) {
    standardGeneric("plot_clusters")
})

setMethod("plot_clusters", "bartsc", function(
    object,
    mod,
    font_size = 11,
    cut = TRUE,
    cut_height,
    gap_size = 5) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (cut == FALSE) { # when it is set to uncut, plot uncut
        return(plot_clusters_uncut(object = object, mod = mod, font_size = font_size))
    } else if (missing(cut_height)) {
        if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_height"]])) { # when no valid cut_height exists, plot uncut
            return(plot_clusters_uncut(object = object, mod = mod, font_size = font_size))
        } else {
            return(plot_clusters_cut(
                object = object, mod = mod, font_size = font_size,
                cut_height = object@resultsCrossCellType[[mod]][["h_clustering"]][["cut_height"]],
                gap_size = gap_size
            ))
        }
    } else {
        return(plot_clusters_cut(object = object, mod = mod, font_size = font_size, cut_height = cut_height, gap_size = gap_size))
    }
})


setGeneric("plot_clusters_uncut", function(
    object,
    mod,
    font_size = 11) {
    standardGeneric("plot_clusters_uncut")
})

#' @importFrom ggdendro dendro_data
#' @importFrom ggdendro segment
#' @importFrom cowplot plot_grid
setMethod("plot_clusters_uncut", "bartsc", function(
    object,
    mod,
    font_size = 11) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]])) {
        stop("No hierarchical clustering result. Run h_clustering(object, mod) first.")
    } else {
        model <- object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]]
    }

    df <- object@resultsCrossCellType[[mod]][["h_clustering"]][["deviation_df"]]
    cell_types_used <- get_used_cell_types(object)
    dims <- length(cell_types_used)

    # dendrogram
    dhc <- as.dendrogram(model)

    plot_data <- dendro_data(dhc, type = "rectangle")
    plot_df1 <- segment(plot_data)
    x_lower <- min(c(plot_df1$x, plot_df1$xend)) - 0.5
    x_upper <- max(c(plot_df1$x, plot_df1$xend)) + 0.5
    y_lower <- min(c(plot_df1$y, plot_df1$yend))
    y_upper <- max(c(plot_df1$y, plot_df1$yend))

    p1 <- ggplot(plot_df1) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_x_continuous(limits = c(x_lower, x_upper), expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_bw() +
        theme(
            text = element_text(size = font_size),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line.x = element_blank(), axis.text.x = element_text(hjust = 1),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            panel.border = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 1), "cm")
        )

    # heatmap
    df_ordered <- df[model$order, ]
    plot_df2 <- expand.grid(X = colnames(df_ordered), Y = rownames(df_ordered))
    values <- c()
    for (i in 1:nrow(plot_df2)) {
        tf <- as.character(plot_df2[i, "Y"])
        comp <- as.character(plot_df2[i, "X"])
        values <- c(values, df_ordered[tf, comp])
    }
    plot_df2$V_ratio <- values

    p2 <- ggplot(plot_df2, aes(X, Y, fill = V_ratio)) +
        geom_tile(width = 1, height = 1) +
        scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = 0, name = "Deviation Ratio") +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme_bw() +
        theme(
            text = element_text(size = font_size),
            panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
        )

    p_out <- plot_grid(p1, p2, align = "h")

    return(p_out)
})

setGeneric("plot_clusters_cut", function(
    object,
    mod,
    font_size = 11,
    cut_height,
    gap_size = 5) {
    standardGeneric("plot_clusters_cut")
})

#' @importFrom ggdendro dendro_data
#' @importFrom ggdendro segment
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
setMethod("plot_clusters_cut", "bartsc", function(
    object,
    mod,
    font_size = 11,
    cut_height,
    gap_size = 5) {
    if (missing(mod)) {
        mod <- get_active_mod(object)
    }

    if (!mod %in% c("RNA", "ATAC", "bimodal")) {
        stop("valid modalities are RNA, ATAC and bimodal")
    }

    if (is.null(object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]])) {
        stop("No hierarchical clustering result. Run h_clustering(object, mod) first.")
    } else {
        model <- object@resultsCrossCellType[[mod]][["h_clustering"]][["model"]]
    }

    df <- object@resultsCrossCellType[[mod]][["h_clustering"]][["deviation_df"]]
    cell_types_used <- get_used_cell_types(object)
    dims <- length(cell_types_used)

    # dendrogram
    dhc <- as.dendrogram(model)

    plot_data <- dendro_data(dhc, type = "rectangle")
    plot_df1 <- segment(plot_data)
    x_lower <- min(c(plot_df1$x, plot_df1$xend)) - 0.5
    x_upper <- max(c(plot_df1$x, plot_df1$xend)) + 0.5
    y_lower <- min(c(plot_df1$y, plot_df1$yend))
    y_upper <- max(c(plot_df1$y, plot_df1$yend))

    p1.1 <- ggplot(plot_df1) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_x_continuous(limits = c(x_lower, x_upper), expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        geom_hline(yintercept = cut_height, linetype = "dashed", color = "blue") +
        theme_bw() +
        theme(
            text = element_text(size = font_size),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line.x = element_blank(), axis.text.x = element_text(hjust = 1),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            panel.border = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 1), "cm")
        )

    cut_clusters <- cutree(model, h = cut_height)

    # rename clusters
    df_ordered <- df[model$order, ]
    cut_clusters <- cut_clusters[rownames(df_ordered)]
    cut_clusters <- factor(cut_clusters, levels = unique(cut_clusters))
    levels(cut_clusters) <- length(levels(cut_clusters)):1
    names(cut_clusters) <- rownames(df_ordered)

    # add blank rows as gaps
    cluster_df <- data.frame(Y = rownames(df_ordered))
    cluster_df$cluster <- cut_clusters[cluster_df$Y]

    ordered_clusters <- unique(cluster_df$cluster) # clusters by the order in heatmap

    cluster_size <- sapply(ordered_clusters, function(x) {
        length(which(cluster_df$cluster == x))
    }) # size of each cluster

    for (i in 1:(length(ordered_clusters) - 1)) {
        head_sec <- ordered_clusters[1:i]
        tail_sec <- ordered_clusters[(i + 1):length(ordered_clusters)]

        head_df <- cluster_df[1:max(which(cluster_df$cluster %in% head_sec)), ]
        tail_df <- cluster_df[min(which(cluster_df$cluster %in% tail_sec)):nrow(cluster_df), ]

        if(gap_size > 0){
            for (j in 1:gap_size) {
                head_df <- rbind(head_df, c(paste0("gap_", i), 0))
            }
        }
        cluster_df <- rbind(head_df, tail_df)
    }

    plot_df2.1 <- expand.grid(X = colnames(df_ordered), Y = cluster_df$Y)
    values <- c()
    for (i in 1:nrow(plot_df2.1)) {
        tf <- as.character(plot_df2.1[i, "Y"])
        comp <- as.character(plot_df2.1[i, "X"])
        if (grepl("gap_\\d+", tf) == TRUE) {
            values <- c(values, NA)
        } else {
            values <- c(values, df_ordered[tf, comp])
        }
    }
    plot_df2.1$V_ratio <- values
    plot_df2.1$y_coord <- rep(seq(1, nrow(cluster_df)), each = dims * (dims - 1))

    lower_y <- min(plot_df2.1$y_coord) - 0.5
    upper_y <- max(plot_df2.1$y_coord) + 0.5

    p2.1 <- ggplot(plot_df2.1[which(!grepl("gap_\\d+", plot_df2.1$Y)), ], aes(X, y_coord, fill = V_ratio)) +
        geom_tile() +
        # scale_fill_gradientn(colors = colorRampPalette(brewer.pal(n = 7, name = "Oranges")) +
        scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = 0, name = "Deviation Ratio") +
        # scale_fill_gradient2(low = "royalblue4", high = "firebrick3", mid = 0, name = "Deviation Ratio") +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        scale_y_continuous(limits = c(lower_y, upper_y), expand = c(0, 0)) +
        theme_bw() +
        theme(
            text = element_text(size = font_size),
            panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
        )

    # update the dendrogram
    plot_df1.1 <- segment(plot_data)

    if(gap_size > 0){
        for (i in 1:(length(ordered_clusters) - 1)) {
            cut_line <- sum(cluster_size[1:i]) + (i - 1) * gap_size + 0.5
            plot_df1.1$x[which(plot_df1.1$x >= cut_line)] <- plot_df1.1$x[which(plot_df1.1$x >= cut_line)] + gap_size
            plot_df1.1$xend[which(plot_df1.1$xend >= cut_line)] <- plot_df1.1$xend[which(plot_df1.1$xend >= cut_line)] + gap_size
        }
    }

    x_lower <- min(c(plot_df1.1$x, plot_df1.1$xend)) - 0.5
    x_upper <- max(c(plot_df1.1$x, plot_df1.1$xend)) + 0.5
    y_lower <- min(c(plot_df1.1$y, plot_df1.1$yend))
    y_upper <- max(c(plot_df1.1$y, plot_df1.1$yend))

    p1.1 <- ggplot(plot_df1.1) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_x_continuous(limits = c(x_lower, x_upper), expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        geom_hline(yintercept = cut_height, linetype = "dashed", color = "blue") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line.x = element_blank(), axis.text.x = element_text(hjust = 1),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            panel.border = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 1), "cm")
        )

    # annotation bar for clusters
    plot_df3.1 <- data.frame(Y = cluster_df$Y, y_coord = 1:nrow(cluster_df), X = 1, cluster = cluster_df$cluster)
    plot_df3.1$cluster <- factor(plot_df3.1$cluster, levels = 1:length(levels(cluster_df$cluster)))

    legend_col <- 1
    if (length(levels(cluster_df$cluster)) > 12) {
        legend_col <- 2
    }

    p3.1 <- ggplot(plot_df3.1[which(!grepl("gap_\\d+", plot_df3.1$Y)), ], aes(X, y_coord, fill = cluster)) +
        geom_tile() +
        scale_y_continuous(limits = c(lower_y, upper_y), expand = c(0, 0)) +
        theme_bw() +
        guides(fill = guide_legend(ncol = legend_col)) +
        theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
        )

    prow <- plot_grid(
        p1.1 + theme(legend.position = "none"),
        p3.1 + theme(legend.position = "none"),
        p2.1 + theme(legend.position = "none"),
        rel_widths = c(20, 1, 20),
        align = "h",
        nrow = 1
    )

    p_out <- plot_grid(
        prow,
        plot_grid(
            get_legend(p2.1),
            get_legend(p3.1),
            NULL,
            ncol = 1,
            align = "v",
            rel_heights = c(1, 2, 1)
        ),
        rel_widths = c(5, 1)
    )

    return(p_out)
})
