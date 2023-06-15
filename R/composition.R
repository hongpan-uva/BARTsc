#' run BART geneset mode
#'
#' This function run BART geneset mode. Input could be a pre-created Bart
#' object, or be set by specifying name, genome and symbol_list
#'
#' @param object a BART object with valid geneset data
#' @param name name of the input data, it could be a sample name, a cluster id and so on
#' @param genome "mm10" or "hg38"
#' @param symbol_list a vector of gene symbols, check the internal data B_cell_gene as an example
#' @param gene_mode_param a list of costumized arguments for gene mode
#'
#' @return A Bart object
#'
#' @export
#'
#' @examples
#' bart_obj <- Bart("B_cell_gene", "hg38", gene_data = B_cell_gene)
#' bart_obj <- runBartGeneSet(bart_obj)
#'
#' bart_obj <- runBartGeneSet("B_cell_gene", "hg38", B_cell_gene)
runBartGeneSet <- function(object, name, genome, symbol_list = NULL,
                           gene_mode_param = list(binsize = 1000)) {
    if (missing(object)) {
        object <- Bart(name, genome, gene_data = symbol_list)
    }
    object <- runMarge(object)
    object <- identifyEnhancer(object)
    object <- predictTF(object, mode = "geneset")
    return(object)
}

#' run BART region mode
#'
#' This function run BART region mode. Input could be a pre-created Bart
#' object, or be set by specifying name, genome and region_list
#'
#' @param object a BART object with valid region data
#' @param name name of the input data, it could be a sample name, a cluster id and so on
#' @param genome "mm10" or "hg38"
#' @param region_list a dataframe that follows a BED6 format, check the internal data B_cell_region as an example
#' @param region_mode_param a list of costumized arguments for region mode
#'
#' @return A Bart object
#'
#' @export
#'
#' @examples
#' bart_obj <- Bart("B_cell_region", "hg38", region_data = B_cell_region)
#' bart_obj <- runBartBartRegion(bart_obj)
#'
#' bart_obj <- runBartBartRegion("B_cell_region", "hg38", B_cell_region)
runBartRegion <- function(object, name, genome, region_list = NULL,
                          gene_mode_param = list(binsize = 1000)) {
    if (missing(object)) {
        object <- Bart(name, genome, region_data = region_list)
    }
    object <- mapRegionScore(object)
    object <- predictTF(object, mode = "region")
    return(object)
}
