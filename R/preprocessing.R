#' Check the existence of a package
#'
#' @param ... Package names
#' @param error If true, throw an error if the package doesn't exist
#'
#' @return Invisibly returns boolean denoting if the package is installed
#'
#' @export
#'
#' @concept utils
#'
#' @section Lifecycle:
#'
#' @examples
#' PackageCheck("SeuratObject", error = FALSE)
#'
PackageCheck <- function(..., error = TRUE) {
    pkgs <- unlist(x = c(...), use.names = FALSE)
    package.installed <- vapply(
        X = pkgs,
        FUN = requireNamespace,
        FUN.VALUE = logical(length = 1L),
        quietly = TRUE
    )
    if (error && any(!package.installed)) {
        stop(
            "Cannot find the following packages: ",
            paste(pkgs[!package.installed], collapse = ", "),
            ". Please install"
        )
    }
    invisible(x = package.installed)
}

WilcoxDETest.old <- function(
    data.use,
    cells.1,
    cells.2,
    verbose = TRUE) {
    data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
    j <- seq_len(length.out = length(x = cells.1))
    overflow.check <- ifelse(
        test = is.na(x = suppressWarnings(length(x = data.use[1, ]) * length(x = data.use[1, ]))),
        yes = FALSE,
        no = TRUE
    )
    limma.check <- PackageCheck("limma", error = FALSE)
    if (limma.check[1] && overflow.check) {
        p_val <- mclapply(
            X = 1:nrow(x = data.use),
            FUN = function(x) {
                return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
            }
        )
    } else {
        if (getOption("BARTsc.limma.wilcox.msg", TRUE) && overflow.check) {
            message(
                "For a more efficient implementation of the Wilcoxon Rank Sum Test,",
                "\n(default method for FindMarkers) please install the limma package",
                "\n--------------------------------------------",
                "\ninstall.packages('BiocManager')",
                "\nBiocManager::install('limma')",
                "\n--------------------------------------------",
                "\nAfter installation of limma, we will automatically use the more ",
                "\nefficient implementation (no further action necessary).",
                "\nThis message will be shown once per session"
            )
            options(BARTsc.limma.wilcox.msg = FALSE)
        }
        group.info <- data.frame(row.names = c(cells.1, cells.2))
        group.info[cells.1, "group"] <- "Group1"
        group.info[cells.2, "group"] <- "Group2"
        group.info[, "group"] <- factor(x = group.info[, "group"])
        data.use <- data.use[, rownames(x = group.info), drop = FALSE]
        p_val <- pblapply(
            X = 1:nrow(x = data.use),
            FUN = function(x) {
                return(wilcox.test(data.use[x, ] ~ group.info[, "group"])$p.value)
            }
        )
    }
    return(data.frame(p_val = unlist(p_val), row.names = rownames(x = data.use)))
}

#' @importFrom presto wilcoxauc
WilcoxDETest <- function(
    data.use,
    cells.1,
    cells.2) {
    label <- colnames(data.use)
    label[which(label %in% cells.1)] <- "fore"
    label[which(label %in% cells.2)] <- "back"
    wilcox_df <- presto::wilcoxauc(data.use, groups_use = c("fore", "back"), label)
    wilcox_df <- wilcox_df[which(wilcox_df$group == "fore"), ]
    return(wilcox_df)
}

# Perform differential expression testing using a logistic regression framework
#
# Constructs a logistic regression model predicting group membership based on a
# given feature and compares this to a null model with a likelihood ratio test.
#
# @param data.use expression matrix
# @param cells.1 a vector of group 1 cell names
# @param cells.2 a vector of group 2 cell names
# @param latent.vars Latent variables to include in model
# @param verbose Print messages
#
#' @importFrom lmtest lrtest
#' @importFrom stats as.formula glm
#' @import pbapply
LRDETest <- function(
    data.use,
    cells.1,
    cells.2,
    latent.vars = NULL,
    verbose = TRUE) {
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(group.info), drop = FALSE]
    latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
    p_val <- pblapply(
        X = 1:nrow(x = data.use),
        FUN = function(x) {
            if (is.null(x = latent.vars)) {
                model.data <- cbind(GENE = data.use[x, ], group.info)
                fmla <- as.formula(object = "group ~ GENE")
                fmla2 <- as.formula(object = "group ~ 1")
            } else {
                model.data <- cbind(GENE = data.use[x, ], group.info, latent.vars)
                fmla <- as.formula(object = paste(
                    "group ~ GENE +",
                    paste(colnames(x = latent.vars), collapse = "+")
                ))
                fmla2 <- as.formula(object = paste(
                    "group ~",
                    paste(colnames(x = latent.vars), collapse = "+")
                ))
            }
            model1 <- glm(formula = fmla, data = model.data, family = "binomial")
            model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
            lrtest <- lrtest(model1, model2)
            return(lrtest$Pr[2])
        }
    )
    to.return <- data.frame(p_val = unlist(p_val), row.names = rownames(data.use))
    return(to.return)
}

#' normalize scATAC-seq count matrix with TF_IDF methods
#'
#' @param method Which TF-IDF implementation to use. Choice of:
#' \itemize{
#'  \item{1}: The TF-IDF implementation used by Stuart & Butler et al. 2019
#'  (\doi{10.1101/460147}). This computes
#'  \eqn{\log(TF \times IDF)}.
#'  \item{2}: The TF-IDF implementation used by Cusanovich & Hill
#'  et al. 2018 (\doi{10.1016/j.cell.2018.06.052}). This
#'  computes \eqn{TF \times (\log(IDF))}.
#'  \item{3}: The log-TF method used by Andrew Hill.
#'  This computes \eqn{\log(TF) \times \log(IDF)}.
#'  \item{4}: The 10x Genomics method (no TF normalization). This computes
#'  \eqn{IDF}.
#' }
#' @param scale.factor Which scale factor to use. Default is 10000.
#' @param idf A precomputed IDF vector to use. If NULL, compute based on the
#' input data matrix.
#' @param verbose Print progress
#' @importFrom Matrix colSums rowSums Diagonal tcrossprod
#' @importFrom methods is "slot<-" slot
#' @export
#' @concept preprocessing
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' RunTFIDF(object = mat)
RunTFIDF <- function(
    object,
    assay = NULL,
    method = 1,
    scale.factor = 1e4,
    idf = NULL,
    verbose = TRUE,
    ...) {
    if (inherits(x = object, what = "data.frame")) {
        object <- as.matrix(x = object)
    }
    if (!inherits(x = object, what = "CsparseMatrix")) {
        object <- as(object = object, Class = "CsparseMatrix")
    }
    if (verbose) {
        message("Performing TF-IDF normalization")
    }
    npeaks <- colSums(x = object)
    if (any(npeaks == 0)) {
        warning("Some cells contain 0 total counts")
    }
    if (method == 4) {
        tf <- object
    } else {
        tf <- tcrossprod(x = object, y = Diagonal(x = 1 / npeaks))
    }
    if (!is.null(x = idf)) {
        precomputed_idf <- TRUE
        if (!inherits(x = idf, what = "numeric")) {
            stop("idf parameter must be a numeric vector")
        }
        if (length(x = idf) != nrow(x = object)) {
            stop(
                "Length of supplied IDF vector does not match",
                " number of rows in input matrix"
            )
        }
        if (any(idf == 0)) {
            stop("Supplied IDF values cannot be zero")
        }
        if (verbose) {
            message("Using precomputed IDF vector")
        }
    } else {
        precomputed_idf <- FALSE
        rsums <- rowSums(x = object)
        if (any(rsums == 0)) {
            warning("Some features contain 0 total counts")
        }
        idf <- ncol(x = object) / rsums
    }

    if (method == 2) {
        if (!precomputed_idf) {
            idf <- log(1 + idf)
        }
    } else if (method == 3) {
        slot(object = tf, name = "x") <- log1p(
            x = slot(object = tf, name = "x") * scale.factor
        )
        if (!precomputed_idf) {
            idf <- log(1 + idf)
        }
    }
    norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
    if (method == 1) {
        slot(object = norm.data, name = "x") <- log1p(
            x = slot(object = norm.data, name = "x") * scale.factor
        )
    }
    colnames(x = norm.data) <- colnames(x = object)
    rownames(x = norm.data) <- rownames(x = object)
    # set NA values to 0
    vals <- slot(object = norm.data, name = "x")
    vals[is.na(x = vals)] <- 0
    slot(object = norm.data, name = "x") <- vals
    return(norm.data)
}
