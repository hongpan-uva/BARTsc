ADAPTIVE_LASSO_MAXSAMPLES <- 20

setClass("Bart",
    slots = c(
        meta = "list",
        param = "list",
        data = "list",
        Args = "list",
        intermediate = "list",
        result = "list"
    ),
    prototype = list(
        param = list(),
        data = list(),
        Args = list(),
        intermediate = list(),
        result = list()
    )
)

Bart <- function(name, genome, gene_data = NULL, region_data = NULL,
                 gene_mode_param = list(binsize = 1000),
                 region_mode_param = list(
                     binsize = 50, scorecol = 5
                 )) {
    # to do: separate generateArgs from helper as addArgs()
    bart_obj <- new("Bart")
    bart_obj@meta$name <- name
    bart_obj@meta$genome <- genome
    # set gene mode data and paramter
    if (!is.null(gene_data)) {
        bart_obj@data[["input_genes"]] <- gene_data
        bart_obj@param[["gene_mode_param"]] <- gene_mode_param
        bart_obj@Args[["gene_mode_args"]] <- generateArgs(bart_obj, "geneset")
    }
    # set region mode data and paramter
    if (!is.null(region_data)) {
        bart_obj@data[["input_regions"]] <- region_data
        bart_obj@param[["region_mode_param"]] <- region_mode_param
        bart_obj@Args[["region_mode_args"]] <- generateArgs(bart_obj, "region")
    }

    return(bart_obj)
}

setGeneric("generateArgs", function(x, subcommand) standardGeneric("generateArgs"))

setMethod("generateArgs", "Bart", function(x, subcommand) {
    # to do: add the bimodal option
    if (subcommand == "geneset") {
        args <- types$SimpleNamespace(
            subcommand_name = "geneset",
            ofilename = x@meta$name,
            species = x@meta$genome,
            binsize = as.integer(x@param[["gene_mode_param"]][["binsize"]]),
            genelist = x@data[["input_genes"]],
            refseq = FALSE, # ignore
            target = NULL, # ignore
            nonorm = FALSE, # ignore
            outdir = "." # ignore
        )
        args <- OptValidator$opt_validate(args)
    } else if (subcommand == "region") {
        args <- types$SimpleNamespace(
            subcommand_name = "region",
            ofilename = x@meta$name,
            species = x@meta$genome,
            binsize = as.integer(x@param[["region_mode_param"]][["binsize"]]),
            scorecol = as.integer(x@param[["region_mode_param"]][["scorecol"]]),
            in_df = x@data[["input_regions"]],
            refseq = FALSE, # ignore
            target = NULL, # ignore
            nonorm = FALSE, # ignore
            outdir = "." # ignore
        )
        args <- OptValidator$opt_validate(args)
    }
    return(args)
})

setGeneric("runMarge", function(x, samples) standardGeneric("runMarge"))

setMethod("runMarge", "Bart", function(x, samples) {
    args <- x@Args[["gene_mode_args"]]
    rp_args <- types$SimpleNamespace(
        genome = args$species,
        histRP = args$rp,
        genelist = args$genelist,
        sym = args$tss,
        name = args$ofilename,
        maxsamples = ADAPTIVE_LASSO_MAXSAMPLES,
        transform = "sqrt",
        exptype = "Gene_Only",
        annotation = args$desc
    )
    marge_res <- RPRegress$main(rp_args, samples)
    names(marge_res) <- c("H3K27ac_selected", "coef")

    x@intermediate[["Marge_based"]][["Marge_model"]] <- marge_res
    return(x)
})

setGeneric("identifyEnhancer", function(x) standardGeneric("identifyEnhancer"))

setMethod("identifyEnhancer", "Bart", function(x) {
    args <- x@Args[["gene_mode_args"]]
    enhancer_args <- types$SimpleNamespace(
        sample_df = x@intermediate[["Marge_based"]][["Marge_model"]][["H3K27ac_selected"]],
        name = args$ofilename,
        k27achdf5 = args$rpkm
    )

    EI_res <- EnhancerIdentifier$main(enhancer_args)

    # get sorted UDHS ids
    positions <- EI_res[[1]] %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    # predicted_enhancers is the 'counting' variable in bart2
    x@intermediate[["Marge_based"]][["predicted_enhancers"]] <- EI_res[[1]]
    x@intermediate[["Marge_based"]][["predicted_enhancers_rank"]] <- positions

    # predicted_enhancers_df is for user to check
    # x@intermediate[["Marge_based"]][["predicted_enhancers_df"]] <- EI_res[[2]]
    return(x)
})

setGeneric("mapRegionScore", function(x) standardGeneric("mapRegionScore"))

setMethod("mapRegionScore", "Bart", function(x) {
    args <- x@Args[["region_mode_args"]]
    counting <- score_on_UDHS$score_on_DHS(args)
    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    x@intermediate[["region_based"]][["overlapped_enhancers"]] <- counting
    x@intermediate[["region_based"]][["overlapped_enhancers_rank"]] <- positions

    return(x)
})

setGeneric("combineModalities", function(x) standardGeneric("combineModalities"))

setMethod("combineModalities", "Bart", function(x) {
    gene_score_li <- x@intermediate[["Marge_based"]][["predicted_enhancers"]]
    gene_score_df <- data.frame(
        ID = as.integer(names(gene_score_li)),
        gene_score = as.numeric(unname(unlist(gene_score_li)))
    ) %>%
        dplyr::arrange(., ID)

    region_score_li <- x@intermediate[["region_based"]][["overlapped_enhancers"]]
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

    x@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
    x@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

    return(x)
})


setGeneric(
    "predictTF",
    function(x, mode) standardGeneric("predictTF")
)

setMethod("predictTF", "Bart", function(x, mode) {
    if (mode == "geneset") {
        INT <- "Marge_based"
        counting <- "predicted_enhancers"
        positions <- "predicted_enhancers_rank"
        args <- x@Args[["gene_mode_args"]]
    } else if (mode == "region") {
        INT <- "region_based"
        counting <- "overlapped_enhancers"
        positions <- "overlapped_enhancers_rank"
        args <- x@Args[["region_mode_args"]]
    } else if (mode == "bimodal") {
        INT <- "bimodal"
        counting <- "combined_enhancers"
        positions <- "combined_enhancers_rank"
        args <- x@Args[["gene_mode_args"]] # temporarily share args with gene mode
    } else {
        stop("wrong mode")
    }
    # find tied intervals
    tied_list <-
        main$find_tied_intervals(
            x@intermediate[[INT]][[counting]],
            x@intermediate[[INT]][[positions]]
        )

    # calculate auc
    auc <-
        AUCcalc$cal_auc(args, x@intermediate[[INT]][[positions]], tied_list)
    auc_df <- data.frame(ChIP_seq = unlist(auc[[2]]), AUC = unlist(auc[[1]]))

    # final stats
    stat_res <- StatTest$stat_test(auc[[1]], auc[[2]], args$normfile)
    tf_names <- attr(stat_res, "row.names")
    stat_res_df <- stat_res %>%
        lapply(., unlist) %>%
        as.data.frame()
    stat_res_df$TR <- tf_names
    stat_res_df <- stat_res_df[, c(
        "TR", "score", "pvalue", "zscore",
        "max_auc", "rank_avg_z_p_a", "rank_avg_z_p_a_irwinhall_pvalue"
    )]

    x@intermediate[[INT]][["tied_list"]] <- tied_list
    x@result[[mode]][["auc"]] <- auc_df
    x@result[[mode]][["stats"]] <- stat_res_df

    return(x)
})
