test_that(".aggregate_rank normalizes equal weights to identical results", {
    rna <- setNames(c(5, 4, 3, 2, 1), as.character(1:5))
    atac <- setNames(c(1, 2, 3, 4, 5), as.character(1:5))

    res1 <- BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(1, 1))
    res2 <- BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(2, 2))
    expect_equal(as.numeric(res1), as.numeric(res2), tolerance = 1e-10)
})

test_that(".aggregate_rank weighted geometric mean is correct", {
    rna <- setNames(c(5, 4, 3, 2, 1), as.character(1:5))
    atac <- setNames(c(1, 2, 3, 4, 5), as.character(1:5))
    w <- c(3, 1)
    w_norm <- w / sum(w)

    res <- BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = w)
    # ranks in UDHS-id order after capping at n_valid=3:
    # RNA: 1, 2, 3, 4.5, 4.5; ATAC: 4.5, 4.5, 3, 2, 1
    expected <- c(
        (1^w_norm[1]) * (4.5^w_norm[2]),
        (2^w_norm[1]) * (4.5^w_norm[2]),
        (3^w_norm[1]) * (3^w_norm[2]),
        (4.5^w_norm[1]) * (2^w_norm[2]),
        (4.5^w_norm[1]) * (1^w_norm[2])
    )
    expect_equal(as.numeric(res), expected, tolerance = 1e-10)
})

test_that(".aggregate_rank weighted MRR is correct", {
    rna <- setNames(c(5, 4, 3, 2, 1), as.character(1:5))
    atac <- setNames(c(1, 2, 3, 4, 5), as.character(1:5))
    w <- c(3, 1)
    w_norm <- w / sum(w)

    res <- BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "MRR", weights = w)
    expected <- c(
        1 / (w_norm[1] / 1 + w_norm[2] / 4.5),
        1 / (w_norm[1] / 2 + w_norm[2] / 4.5),
        1 / (w_norm[1] / 3 + w_norm[2] / 3),
        1 / (w_norm[1] / 4.5 + w_norm[2] / 2),
        1 / (w_norm[1] / 4.5 + w_norm[2] / 1)
    )
    expect_equal(as.numeric(res), expected, tolerance = 1e-10)
})

test_that(".aggregate_rank weighted arithmetic mean is correct", {
    rna <- setNames(c(5, 4, 3, 2, 1), as.character(1:5))
    atac <- setNames(c(1, 2, 3, 4, 5), as.character(1:5))
    w <- c(3, 1)
    w_norm <- w / sum(w)

    res <- BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "mean", weights = w)
    expected <- c(
        w_norm[1] * 1 + w_norm[2] * 4.5,
        w_norm[1] * 2 + w_norm[2] * 4.5,
        w_norm[1] * 3 + w_norm[2] * 3,
        w_norm[1] * 4.5 + w_norm[2] * 2,
        w_norm[1] * 4.5 + w_norm[2] * 1
    )
    expect_equal(as.numeric(res), expected, tolerance = 1e-10)
})

test_that(".aggregate_rank rejects invalid weights", {
    rna <- setNames(c(5, 4, 3, 2, 1), as.character(1:5))
    atac <- setNames(c(1, 2, 3, 4, 5), as.character(1:5))

    expect_error(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(1, -1)))
    expect_error(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(1, 0)))
    expect_error(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(1)))
    expect_error(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(1, 1, 1)))
    expect_error(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c(NA, 1)))
    expect_error(suppressWarnings(BARTsc:::.aggregate_rank(rna, atac, n_valid = 3, method = "geom.mean", weights = c("a", "b"))))
})

test_that("bartsc stores custom bimodal weights", {
    obj <- bartsc(
        name = "test",
        genome = "hg38",
        label = factor(c("A", "B")),
        RNA_cnt_matrix = matrix(0, nrow = 1, ncol = 2),
        ATAC_cnt_matrix = Matrix::Matrix(0, nrow = 1, ncol = 2, sparse = TRUE),
        peaks = data.frame(seqnames = "chr1", start = 1, end = 2, name = ".", score = ".", strand = "."),
        bimodal_mode_param = list(binsize = 50, weights = c(3, 1))
    )
    expect_equal(obj@param$bimodal_mode_param$weights, c(3, 1))
})

test_that("combineModsByTopRank reads weights from bart object", {
    obj <- new("bart")
    obj@meta$name <- "test"
    obj@meta$genome <- "hg38"
    obj@param$bimodal_mode_param <- list(binsize = 50, weights = c(3, 1))
    obj@intermediate$Marge_based$predicted_enhancers <- as.list(setNames(c(5, 4, 3, 2, 1), as.character(1:5)))
    obj@intermediate$region_based$overlapped_enhancers <- as.list(setNames(c(1, 2, 3, 4, 5), as.character(1:5)))

    obj <- BARTsc:::combineModsByTopRank(obj, n_valid = 3, method = "geom.mean")

    w_norm <- c(3, 1) / sum(c(3, 1))
    expected <- c(
        (1^w_norm[1]) * (4.5^w_norm[2]),
        (2^w_norm[1]) * (4.5^w_norm[2]),
        (3^w_norm[1]) * (3^w_norm[2]),
        (4.5^w_norm[1]) * (2^w_norm[2]),
        (4.5^w_norm[1]) * (1^w_norm[2])
    )
    expect_equal(as.numeric(unlist(obj@intermediate$bimodal$combined_enhancers)), -expected, tolerance = 1e-10)
})
