ADAPTIVE_LASSO_MAXSAMPLES <- 20

#' An S4 class to store raw data, intermediate data and outcomes of BARTsc analysis
#'
#' @slot meta, param, data, Args, intermediate, result
setClass("bart",
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

#' helper of bart class
#'
#' This function creates a bart object with specified data.
#'
#' @param name name of the input data, it could be a sample name, a cluster id and so on
#' @param genome "mm10" or "hg38"
#' @param gene_data a vector of gene symbols, check the internal data B_cell_gene as an example
#' @param region_data a dataframe that follows a BED6 format, check the internal data B_cell_region as an example
#' @param gene_mode_param a list of costumized arguments for gene mode
#' @param region_mode_param a list of costumized arguments for region mode
#'
#' @return A bart object
#'
#' @export
#'
bart <- function(name, genome, gene_data, region_data,
                 gene_mode_param = list(binsize = 1000),
                 region_mode_param = list(
                     binsize = 50, scorecol = 5
                 ),
                 bimodal_mode_param = list(
                     binsize = 50
                 )) {
    # to do: separate generateArgs from helper as addArgs()
    bart_obj <- new("bart")
    bart_obj@meta$name <- name
    bart_obj@meta$genome <- genome
    # set gene mode data and paramter
    if (!missing(gene_data)) {
        bart_obj@data <- c(bart_obj@data, list(input_genes = gene_data))
        bart_obj@param[["gene_mode_param"]] <- gene_mode_param
        bart_obj@Args[["gene_mode_args"]] <- generateArgs(bart_obj, "geneset")
    }
    # set region mode data and paramter
    if (!missing(region_data)) {
        bart_obj@data <- c(bart_obj@data, list(input_regions = region_data))
        bart_obj@param[["region_mode_param"]] <- region_mode_param
        bart_obj@Args[["region_mode_args"]] <- generateArgs(bart_obj, "region")
    }
    if (!missing(gene_data) && !missing(region_data)) {
        bart_obj@param[["bimodal_mode_param"]] <- bimodal_mode_param
        bart_obj@Args[["bimodal_mode_args"]] <- generateArgs(bart_obj, "bimodal")
    }

    return(bart_obj)
}

#' Create a Args Namespace
#'
#' This function create a python Namespace that holds BART2 arguments.
#'
#' @param object a bart object
#' @param subcommand "geneset" or "region"
#'
#' @return a python Namespace
setGeneric("generateArgs", function(object, subcommand) standardGeneric("generateArgs"))

setMethod("generateArgs", "bart", function(object, subcommand) {
    # to do: add the bimodal option
    if (subcommand == "geneset") {
        args <- types$SimpleNamespace(
            subcommand_name = "geneset",
            ofilename = object@meta$name,
            species = object@meta$genome,
            binsize = as.integer(object@param[["gene_mode_param"]][["binsize"]]),
            genelist = object@data[["input_genes"]],
            refseq = FALSE, # ignore
            target = NULL, # ignore
            nonorm = FALSE, # ignore
            outdir = "." # ignore
        )
        args <- OptValidator$opt_validate(args)
    } else if (subcommand == "region") {
        args <- types$SimpleNamespace(
            subcommand_name = "region",
            ofilename = object@meta$name,
            species = object@meta$genome,
            binsize = as.integer(object@param[["region_mode_param"]][["binsize"]]),
            scorecol = as.integer(object@param[["region_mode_param"]][["scorecol"]]),
            in_df = object@data[["input_regions"]],
            refseq = FALSE, # ignore
            target = NULL, # ignore
            nonorm = FALSE, # ignore
            outdir = "." # ignore
        )
        args <- OptValidator$opt_validate(args)
    } else if (subcommand == "bimodal") {
        args <- types$SimpleNamespace(
            subcommand_name = "region", # use region mode normfile
            ofilename = object@meta$name,
            species = object@meta$genome,
            binsize = as.integer(object@param[["bimodal_mode_param"]][["binsize"]]),
            refseq = FALSE, # ignore
            target = NULL, # ignore
            nonorm = FALSE, # ignore
            outdir = "." # ignore
        )
        args <- OptValidator$opt_validate(args)
    }
    return(args)
})

#' Run Marge on geneset input
#'
#' This function generate a Marge model based on input gene list
#'
#' @param object a bart object
#'
#' @return a bart object
setGeneric("runMarge", function(object) standardGeneric("runMarge"))

setMethod("runMarge", "bart", function(object) {
    args <- object@Args[["gene_mode_args"]]
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
    marge_res <- RPRegress$main(rp_args)
    names(marge_res) <- c("H3K27ac_selected", "coef")

    object@intermediate[["Marge_based"]][["Marge_model"]] <- marge_res
    return(object)
})

#' Identify active enhancers
#'
#' This function generate a cis-regulatory profile using the linear combination of H3K27ac profiles predicted by Marge.
#'
#' @param object a bart object
#'
#' @return a bart object
setGeneric("identifyEnhancer", function(object) standardGeneric("identifyEnhancer"))

setMethod("identifyEnhancer", "bart", function(object) {
    args <- object@Args[["gene_mode_args"]]
    enhancer_args <- types$SimpleNamespace(
        sample_df = object@intermediate[["Marge_based"]][["Marge_model"]][["H3K27ac_selected"]],
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
    object@intermediate[["Marge_based"]][["predicted_enhancers"]] <- EI_res[[1]]
    object@intermediate[["Marge_based"]][["predicted_enhancers_rank"]] <- positions

    # predicted_enhancers_df is for user to check
    # object@intermediate[["Marge_based"]][["predicted_enhancers_df"]] <- EI_res[[2]]
    return(object)
})

#' Map region scores to UDHS sites
#'
#' This function Map region scores to UDHS sites to generate a region-derived cis-regulatory profile
#'
#' @param object a bart object
#'
#' @return a bart object
setGeneric("mapRegionScore", function(object) standardGeneric("mapRegionScore"))

setMethod("mapRegionScore", "bart", function(object) {
    args <- object@Args[["region_mode_args"]]
    counting <- score_on_UDHS$score_on_DHS(args)
    # get sorted UDHS ids
    positions <- counting %>%
        unlist() %>%
        sort(decreasing = TRUE) %>%
        names()

    object@intermediate[["region_based"]][["overlapped_enhancers"]] <- counting
    object@intermediate[["region_based"]][["overlapped_enhancers_rank"]] <- positions

    return(object)
})

setGeneric(
    "predictTF",
    function(object, mode, reserve_interm = FALSE) standardGeneric("predictTF")
)

setMethod("predictTF", "bart", function(object, mode, reserve_interm = FALSE) {
    if (mode == "geneset") {
        INT <- "Marge_based"
        counting <- "predicted_enhancers"
        positions <- "predicted_enhancers_rank"
        args <- object@Args[["gene_mode_args"]]
    } else if (mode == "region") {
        INT <- "region_based"
        counting <- "overlapped_enhancers"
        positions <- "overlapped_enhancers_rank"
        args <- object@Args[["region_mode_args"]]
    } else if (mode == "bimodal") {
        INT <- "bimodal"
        counting <- "combined_enhancers"
        positions <- "combined_enhancers_rank"
        args <- object@Args[["bimodal_mode_args"]] # temporarily share args with gene mode
    } else {
        stop("wrong mode")
    }
    # find tied intervals
    tied_list <-
        main$find_tied_intervals(
            object@intermediate[[INT]][[counting]],
            object@intermediate[[INT]][[positions]]
        )

    print(paste(length(tied_list), "tied intervals"))

    # calculate auc
    auc <-
        AUCcalc$cal_auc(args, object@intermediate[[INT]][[positions]], tied_list)
    auc_df <- data.frame(ChIP_seq = unlist(auc[[2]]), AUC = unlist(auc[[1]]))

    # final stats
    stat_res <- StatTest$stat_test(auc[[1]], auc[[2]], args$normfile)
    tf_names <- attr(stat_res, "row.names")
    stat_res_df <- stat_res %>%
        lapply(., unlist) %>%
        as.data.frame()
    stat_res_df$TF <- tf_names
    stat_res_df <- stat_res_df[, c(
        "TF", "avg_auc", "score", "pvalue", "zscore",
        "max_auc", "rank_avg_z_p_a", "rank_avg_z_p_a_irwinhall_pvalue"
    )]

    if (reserve_interm == TRUE) {
        object@intermediate[[INT]][["tied_list"]] <- tied_list
    } else {
        object@intermediate <- list()
    }
    object@result[[mode]][["auc"]] <- auc_df
    object@result[[mode]][["stats"]] <- stat_res_df

    return(object)
})


setGeneric(
    "generate_null",
    function(object, mode, reserve_interm = FALSE) standardGeneric("generate_null")
)

setMethod("generate_null", "bart", function(object, mode, reserve_interm = FALSE) {
    if (mode == "geneset") {
        INT <- "Marge_based"
        counting <- "predicted_enhancers"
        positions <- "predicted_enhancers_rank"
        args <- object@Args[["gene_mode_args"]]
    } else if (mode == "region") {
        INT <- "region_based"
        counting <- "overlapped_enhancers"
        positions <- "overlapped_enhancers_rank"
        args <- object@Args[["region_mode_args"]]
    } else if (mode == "bimodal") {
        INT <- "bimodal"
        counting <- "combined_enhancers"
        positions <- "combined_enhancers_rank"
        args <- object@Args[["bimodal_mode_args"]] # temporarily share args with gene mode
    } else {
        stop("wrong mode")
    }

    utils::data(null_auc)
    auc_df <- data.frame(ChIP_seq = unlist(auc[[2]]), AUC = unlist(auc[[1]]))

    # final stats
    stat_res <- StatTest$stat_test(auc[[1]], auc[[2]], args$normfile)
    tf_names <- attr(stat_res, "row.names")
    stat_res_df <- stat_res %>%
        lapply(., unlist) %>%
        as.data.frame()
    stat_res_df$TF <- tf_names
    stat_res_df <- stat_res_df[, c(
        "TF", "avg_auc", "score", "pvalue", "zscore",
        "max_auc", "rank_avg_z_p_a", "rank_avg_z_p_a_irwinhall_pvalue"
    )]

    object@result[[mode]][["auc"]] <- auc_df
    object@result[[mode]][["stats"]] <- stat_res_df

    return(object)
})
