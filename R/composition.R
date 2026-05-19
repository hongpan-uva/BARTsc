#' run BART geneset mode
#'
#' This function run BART geneset mode. Input could be a pre-created bart
#' object, or be set by specifying name, genome and symbol_list
#'
#' @param object a BART object with valid geneset data
#' @param name name of the input data, it could be a sample name, a cluster id
#' and so on
#' @param genome "mm10" or "hg38"
#' @param symbol_list a vector of gene symbols, check the internal data
#' B_cell_gene as an example
#' @param gene_mode_param a list of costumized arguments for gene mode
#'
#' @return A bart object
#'
#' @export
#'
run_bart_gene_set <- function(object, name, genome, symbol_list = NULL,
                              gene_mode_param = list(binsize = 1000),
                              reserve_interm = FALSE, return_null = FALSE) {
    if (missing(object)) {
        object <- bart(name, genome,
            gene_data = symbol_list,
            gene_mode_param = gene_mode_param
        )
    }

    if (return_null == TRUE) {
        object <- generate_null(object, mode = "geneset", reserve_interm)
        return(object)
    }

    object <- runMarge(object)
    object <- identifyEnhancer(object)

    # profile_1 <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])

    # # keep the score of top n and bottom n
    # n_top <- 2000
    # n_bottom <- 0
    # profile_1[(n_top + 1):(length(profile_1) - n_bottom)] <- 0
    # object@intermediate[["Marge_based"]][["predicted_enhancers"]] <- as.list(profile_1)
    # object@intermediate[["Marge_based"]][["predicted_enhancers_rank"]] <- profile_1 %>%
    #     unlist() %>%
    #     sort(decreasing = TRUE) %>%
    #     names()

    object <- predictTF(object, mode = "geneset", reserve_interm)
    return(object)
}

#' run BART region mode
#'
#' This function run BART region mode. Input could be a pre-created bart
#' object, or be set by specifying name, genome and region_list
#'
#' @param object a BART object with valid region data
#' @param name name of the input data, it could be a sample name, a cluster id
#' and so on
#' @param genome "mm10" or "hg38"
#' @param region_list a dataframe that follows a BED6 format, check the
#' internal data B_cell_region as an example
#' @param region_mode_param a list of costumized arguments for region mode
#'
#' @return A bart object
#'
#' @export
#'
run_bart_region <- function(object, name, genome, region_list = NULL,
                            region_mode_param = list(
                                binsize = 50, scorecol = 5
                            ),
                            reserve_interm = FALSE, return_null = FALSE) {
    if (missing(object)) {
        object <- bart(name, genome,
            region_data = region_list,
            region_mode_param = region_mode_param
        )
    }

    if (return_null == TRUE) {
        object <- generate_null(object, mode = "region", reserve_interm)
        return(object)
    }

    object <- mapRegionScore(object)
    object <- predictTF(object, mode = "region", reserve_interm)
    return(object)
}

#' run BART bimodal mode
#'
#' This function run BART bimodal mode. Input could be a pre-created bart
#' object, or be set by specifying name, genome, symbol_list and region_list
#'
#' @param object a BART object with valid geneset data
#' @param name name of the input data, it could be a sample name, a cluster id
#' and so on
#' @param genome "mm10" or "hg38"
#' @param symbol_list a vector of gene symbols, check the internal data
#' B_cell_gene as an example
#' @param region_list a dataframe that follows a BED6 format, check the
#' internal data B_cell_region as an example
#' @param gene_mode_param a list of costumized arguments for gene mode
#'
#' @return A bart object
#'
#' @import mgcv
#' @import gratia
#' @export
#'
run_bart_bimodal <- function(
    object, name, genome,
    symbol_list = NULL, region_list = NULL,
    gene_mode_param = list(binsize = 1000),
    region_mode_param = list(
        binsize = 50, scorecol = 5
    ),
    bimodal_mode_param = list(
        binsize = 50
    ),
    reserve_interm = FALSE, return_null = FALSE) {
    DFLT_INT_NUM <- 1000 # default integration number
    RNA_ONLY <- FALSE # whether only use RNA side
    ATAC_ONLY <- FALSE # whether only use ATAC side

    if (missing(object)) {
        object <- bart(name, genome,
            gene_data = symbol_list, region_data = region_list,
            gene_mode_param = gene_mode_param,
            region_mode_param = region_mode_param,
            bimodal_mode_param = bimodal_mode_param
        )
    }

    if (return_null == TRUE) {
        object <- generate_null(object, mode = "bimodal", reserve_interm)
        return(object)
    }

    # predict cis-regulatory profile
    RNA.res <- tryCatch(
        {
            object <- runMarge(object)
            object <- identifyEnhancer(object)
        },
        error = function(cond) {
            return(NA)
        }
    )

    if (is.na(RNA.res)) {
        ATAC_ONLY <- TRUE
    }

    object <- mapRegionScore(object)

    udhs_ATAC_score <- object@intermediate$region_based$overlapped_enhancers
    n_positive <- length(udhs_ATAC_score[which(udhs_ATAC_score > 0)])

    if (n_positive == 0) {
        RNA_ONLY <- TRUE
    }

    if (RNA_ONLY == TRUE && ATAC_ONLY == TRUE) {
        message("Warning: Both RNA and ATAC did not map to enough UDHSs, return null output")
        object <- run_bart_bimodal(object, return_null = TRUE, reserve_interm = reserve_interm)
    } else if (RNA_ONLY == FALSE && ATAC_ONLY == FALSE) {
        udhs_RNA_ordered <- object@intermediate$Marge_based$predicted_enhancers_rank
        udhs_ATAC_ordered <- object@intermediate$region_based$overlapped_enhancers_rank

        overlap_count <- sapply(seq(0, n_positive, by = 10), function(x) {
            intersect(udhs_RNA_ordered[1:x], udhs_ATAC_ordered[1:x]) %>%
                length() %>%
                return()
        })

        df <- data.frame(top_sites_included = seq(0, n_positive, by = 10), overlap_count)
        colnames(df) <- c("x", "y")
        use_default <- FALSE
        out <- tryCatch(
            {
                gam_model <- mgcv::gam(formula = y ~ s(x, bs = "cs"), data = df, method = "REML")
                deriv_df <- gratia::derivatives(gam_model)
                integ_num <- round(deriv_df$x[which.max(deriv_df$.derivative)])
                if (integ_num >= DFLT_INT_NUM) {
                    message(paste(object@meta$name, "top udhs included:", integ_num))
                } else {
                    message(paste(object@meta$name, "optimized integration number too small:", integ_num))
                    integ_num <- DFLT_INT_NUM
                    message(paste("Use default value", integ_num))
                }
            },
            error = function(cond) {
                message(paste("failed to find optimized integration number for ", object@meta$name))
                message(paste("Use default value", DFLT_INT_NUM))
                return("use default")
            }
        )

        if (!is.null(out)) {
            integ_num <- DFLT_INT_NUM
        }

        object <- combineModsByTopRank(object, n_valid = integ_num, method = "geom.mean")
        # object <- combineModsBySign(object)
        # object <- combineModsByProduct(object, ATAC_weight = 1)
        # object <- refineRNAWithATAC(object)
        object <- predictTF(object, "bimodal", reserve_interm)
        message(paste0("integrated with top ", integ_num, " UDHSs"))
    } else if (RNA_ONLY == TRUE) {
        message("Warning: ATAC input mapped to 0 region, will only use RNA side")

        profile_used <- unlist(object@intermediate[["Marge_based"]][["predicted_enhancers"]])
        profile_used[(DFLT_INT_NUM + 1):length(profile_used)] <- -Inf # only use top DFLT_INT_NUM UDHS
        counting <- as.list(profile_used)

        # get sorted UDHS ids
        positions <- counting %>%
            unlist() %>%
            sort(decreasing = TRUE) %>%
            names()

        object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
        object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

        object <- predictTF(object, "bimodal", reserve_interm)
    } else if (ATAC_ONLY == TRUE) {
        message("Warning: RNA side failed, will only use ATAC side")

        profile_used <- unlist(object@intermediate[["region_based"]][["overlapped_enhancers"]])
        profile_used[(DFLT_INT_NUM + 1):length(profile_used)] <- -Inf # only use top DFLT_INT_NUM UDHS
        counting <- as.list(profile_used)

        # get sorted UDHS ids
        positions <- counting %>%
            unlist() %>%
            sort(decreasing = TRUE) %>%
            names()

        object@intermediate[["bimodal"]][["combined_enhancers"]] <- counting
        object@intermediate[["bimodal"]][["combined_enhancers_rank"]] <- positions

        object <- predictTF(object, "bimodal", reserve_interm)
    }

    return(object)
}

#' Safe wrapper for run_bart_gene_set
#'
#' Runs BART geneset mode and falls back to a null result on empty input
#' or any runtime error (e.g., Python sklearn failures).
#'
#' @param object a BART object with valid geneset data
#' @param reserve_interm whether to reserve intermediate results
#'
#' @return A bart object
#'
#' @noRd
safe_run_bart_gene_set <- function(object, reserve_interm = FALSE) {
    if (length(object@data$input_genes) == 0) {
        return(run_bart_gene_set(object, reserve_interm = reserve_interm, return_null = TRUE))
    }
    out <- tryCatch(
        run_bart_gene_set(object, reserve_interm = reserve_interm),
        error = function(e) {
            warning(paste0("BART geneset failed for '", object@meta$name,
                   "': ", conditionMessage(e), ". Returning null result."))
            run_bart_gene_set(object, reserve_interm = reserve_interm, return_null = TRUE)
        }
    )
    return(out)
}

#' Safe wrapper for run_bart_region
#'
#' Runs BART region mode and falls back to a null result on empty input
#' or any runtime error.
#'
#' @param object a BART object with valid region data
#' @param reserve_interm whether to reserve intermediate results
#'
#' @return A bart object
#'
#' @noRd
safe_run_bart_region <- function(object, reserve_interm = FALSE) {
    if (is.null(object@data$input_regions) || nrow(object@data$input_regions) == 0) {
        return(run_bart_region(object, reserve_interm = reserve_interm, return_null = TRUE))
    }
    out <- tryCatch(
        run_bart_region(object, reserve_interm = reserve_interm),
        error = function(e) {
            warning(paste0("BART region failed for '", object@meta$name,
                   "': ", conditionMessage(e), ". Returning null result."))
            run_bart_region(object, reserve_interm = reserve_interm, return_null = TRUE)
        }
    )
    return(out)
}

#' Safe wrapper for run_bart_bimodal
#'
#' Runs BART bimodal mode and falls back to a null result on empty input
#' or any runtime error.
#'
#' @param object a BART object with valid geneset and/or region data
#' @param reserve_interm whether to reserve intermediate results
#'
#' @return A bart object
#'
#' @noRd
safe_run_bart_bimodal <- function(object, reserve_interm = FALSE) {
    gene_empty <- is.null(object@data$input_genes) || length(object@data$input_genes) == 0
    region_empty <- is.null(object@data$input_regions) || nrow(object@data$input_regions) == 0
    if (gene_empty && region_empty) {
        return(run_bart_bimodal(object, reserve_interm = reserve_interm, return_null = TRUE))
    }
    out <- tryCatch(
        run_bart_bimodal(object, reserve_interm = reserve_interm),
        error = function(e) {
            warning(paste0("BART bimodal failed for '", object@meta$name,
                   "': ", conditionMessage(e), ". Returning null result."))
            run_bart_bimodal(object, reserve_interm = reserve_interm, return_null = TRUE)
        }
    )
    return(out)
}
