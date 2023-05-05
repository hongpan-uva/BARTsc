.select_feature <- function(
    matrix1 = NULL,
    prop = 0.25) {
    cutoff_ncell_1 <- ncol(matrix1) * prop
    # select genes expressed in at least 25% cells in group 1
    select_tag_1 <- apply(matrix1, 1, function(x) {
        return(length(x[x > 0]) >= cutoff_ncell_1)
    })
    return(names(select_tag_1)[which(select_tag_1 == TRUE)])
}

get_single_group_markers <- function(
    data_matrix = NULL,
    meta_matrix = NULL,
    col_use = NULL,
    group_1 = NULL,
    group_2 = NULL,
    test_method = "wilcoxon",
    prop = 0.25,
    pvalue_thres = NULL,
    pvalue_adj_thres = NULL,
    log2fc_thres = NULL) {
    groups <- meta_matrix[, col_use] %>%
        table() %>%
        names()

    if (!group_1 %in% groups) {
        stop("Wrong group 1 name")
    }

    if (missing(group_2)) {
        if (!group_2 %in% groups) {
            stop("Wrong group 2 name")
        }
    } else {
        group_2 <- groups[which(groups != group_1)]
    }

    names_fore <- rownames(meta_matrix)[which(meta_matrix[, col_use] %in% group_1)]
    names_back <- rownames(meta_matrix)[which(meta_matrix[, col_use] %in% group_2)]

    data_matrix_fore <- data_matrix[, names_fore]
    data_matrix_back <- data_matrix[, names_back]

    genes_select_1 <- .select_feature(data_matrix_fore, prop)
    data_matrix_fore_use <- data_matrix_fore[genes_select_1, ]
    data_matrix_back_use <- data_matrix_back[genes_select_1, ]

    log2mean_fore <- apply(data_matrix_fore_use, 1, function(y) {
        return(log2(mean(y) + 1))
    })
    log2mean_back <- apply(data_matrix_back_use, 1, function(y) {
        return(log2(mean(y) + 1))
    })

    log2fc <- log2mean_fore - log2mean_back

    # only consider postively different features
    genes_select_2 <- genes_select_1[which(log2fc > 0)]
    log2fc_select_2 <- log2fc[genes_select_2]

    # calculate p.value
    pvalue <- sapply(genes_select_2, function(y) {
        test_out <- wilcox.test(data_matrix_fore[y, ], data_matrix_back[y, ])
        return(test_out$p.value)
    })

    pvalue_adj <- p.adjust(pvalue,
        method = "bonferroni",
        n = nrow(data_matrix_fore)
    )

    out_df <- data.frame(
        log2fc = log2fc_select_2,
        p.value = pvalue, p.value.adj = pvalue_adj
    ) %>%
        {
            if (!is.null(pvalue_thres)) dplyr::filter(., p.value <= pvalue_thres) else .
        } %>%
        {
            if (!is.null(log2fc_thres)) dplyr::filter(., log2fc >= log2fc_thres) else .
        } %>%
        {
            if (!is.null(pvalue_adj_thres)) {
                dplyr::filter(., p.value.adj <= pvalue_adj_thres)
            } else {
                .
            }
        }

    out_df <- out_df[order(out_df$p.value.adj), ]

    return(out_df)
}
