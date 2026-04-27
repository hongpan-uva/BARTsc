#' try to load bart2 and modules of it
#'
#' @return no return
#'
.try_import <- function() {
    tryCatch(
        {
            bart2 <<- reticulate::import("bart2", delay_load = TRUE)
            message("bart2 has been loaded successfully and is now ready for use!")
        },
        error = function(err) {
            stop("No available bart2", call. = FALSE)
        }
    )

    contents <- reticulate::py_get_attr(bart2, "__all__") %>% reticulate::py_to_r()

    for (module in contents) {
        paste0("bart2.", module) %>%
            reticulate::import(., delay_load = TRUE) %>%
            assign(module, value = ., env = .GlobalEnv)
    }
}

#' load installed bart2 and modules of it
#'
#' @param condaenv conda environemnt name
#' @param virtualenv virtual envrionment name
#'
#' @return no return
#'
#' @export
#'
load_bart2 <- function(
    condaenv = NULL,
    virtualenv = "bartsc_env") {
    # load bart2
    if (!is.null(condaenv)) {
        reticulate::use_condaenv(condaenv)
        .try_import()
    } else {
        reticulate::use_virtualenv(virtualenv)
        .try_import()
    }
    # load other python packages to use
    types <<- reticulate::import("types")
}

#' install bart2 python package
#'
#' @param bart_dir path to install bart2
    #' @param python_ver python version to use, default is "3.9"
    #' @param lib path of bart2 library data
#'
#' @return no return
#'
#' @export
#'
install_bart2 <- function(
    bart_dir = NULL,
    python_ver = "3.9",
    lib = NULL) {
    # Ensure the requested Python version is installed and capture its path
    tryCatch(
        {
            reticulate::use_python_version(python_ver)
        },
        error = function(cond) {
            message("Given python version was not found, will install it")
            reticulate::install_python(python_ver)
            reticulate::use_python_version(python_ver)
        }
    )
    python_path <- reticulate::py_exe()

    # Ensure a clean virtual environment
    tryCatch({
        reticulate::virtualenv_remove("bartsc_env", confirm = FALSE)
    }, error = function(e) invisible(NULL))

    # create a virtual env using the exact Python binary
    reticulate::virtualenv_create(envname = "bartsc_env", python = python_path)
    reticulate::use_virtualenv("bartsc_env")

    if (is.null(bart_dir) || bart_dir == "") {
        bart_dir <- fs::path(.libPaths()[[1]], "BARTsc")
    }

    bart_dir <- fs::path_norm(bart_dir)
    fs::dir_create(bart_dir, recurse = TRUE)

    bart2_python_dir <- fs::path(bart_dir, "bart2_python")

    # remove existing bart2
    if (fs::dir_exists(bart2_python_dir)) {
        fs::dir_delete(bart2_python_dir)
    }

    tgz_path <- fs::path(bart_dir, "bart2_python.tgz")
    tmp_extract <- fs::path(bart_dir, "bart2_extract_tmp")

    # register cleanup in case of failure
    success <- FALSE
    on.exit({
        if (!success) {
            if (fs::file_exists(tgz_path)) fs::file_delete(tgz_path)
            if (fs::dir_exists(tmp_extract)) fs::dir_delete(tmp_extract)
            if (fs::dir_exists(bart2_python_dir)) fs::dir_delete(bart2_python_dir)
        }
    }, add = TRUE)

    # download latest bart2
    tryCatch({
        download.file(
            url = "https://github.com/hongpan-uva/bart2/tarball/use_in_r",
            destfile = tgz_path,
            mode = "wb"
        )
    }, error = function(e) {
        stop("Failed to download bart2 source code: ", conditionMessage(e), call. = FALSE)
    })

    # extract tarball
    if (fs::dir_exists(tmp_extract)) {
        fs::dir_delete(tmp_extract)
    }
    fs::dir_create(tmp_extract, recurse = TRUE)

    tryCatch({
        untar(tgz_path, exdir = tmp_extract)
    }, error = function(e) {
        stop("Failed to extract bart2 archive: ", conditionMessage(e), call. = FALSE)
    })

    # identify extracted directory and rename to bart2_python
    extracted_items <- fs::dir_ls(tmp_extract, type = "directory")
    if (length(extracted_items) == 0) {
        stop("Failed to extract bart2_python archive", call. = FALSE)
    }

    fs::file_move(extracted_items[1], bart2_python_dir)
    fs::dir_delete(tmp_extract)
    fs::file_delete(tgz_path)

    # set up config
    conf_file <- fs::path(bart2_python_dir, "bart2", "bart.conf")
    if (fs::file_exists(conf_file) && !is.null(lib)) {
        lines <- readLines(conf_file)
        lines <- sub("=.*", paste0("= ", lib), lines)
        writeLines(lines, conf_file)
    }

    # install bart2 python package
    reticulate::py_install(bart2_python_dir, pip = TRUE)

    success <- TRUE
}

#' Download bart2 library
#'
#' @param lib_dir path to store library data
#' @param site download site, supports "box" or "zenodo"
#'
#' @return path where the library date is stored
#'
#' @export
#'
get_library <- function(lib_dir, site = "box") {
    # select download site for library data
    if (!site %in% c("box", "zenodo")) {
        stop("Invalid site: must be 'box' or 'zenodo'", call. = FALSE)
    }

    if (lib_dir == "") {
        lib_dir <- fs::path(.libPaths()[[1]], "BARTsc")
    }

    lib_dir <- fs::path_norm(lib_dir)
    fs::dir_create(lib_dir, recurse = TRUE)

    lib_path <- fs::path(lib_dir, "bart2_library")

    if (fs::dir_exists(lib_path)) {
        return(lib_path)
    }

    fs::dir_create(lib_path, recurse = TRUE)

    # download links for different sites
    if (site == "box") {
        hg38_url <- "https://virginia.box.com/shared/static/2kqczz9gixetcr9p4bl650uyrio5zd33.gz"
        mm10_url <- "https://virginia.box.com/shared/static/bxdggnhp4bjz2l5h2zjlisnzp0ac7axf.gz"
    } else {
        # pseudo links for zenodo
        hg38_url <- "https://zenodo.org/records/18854649/files/hg38_library.tar.gz?download=1"
        mm10_url <- "https://zenodo.org/records/18854649/files/mm10_library.tar.gz?download=1"
    }

    # download and extract hg38
    hg38_tgz <- fs::path(lib_path, "hg38_library.tar.gz")
    tryCatch({
        download.file(hg38_url, destfile = hg38_tgz, mode = "wb")
        untar(hg38_tgz, exdir = lib_path)
    }, error = function(e) {
        if (fs::file_exists(hg38_tgz)) fs::file_delete(hg38_tgz)
        stop("Failed to download or extract hg38 library: ", conditionMessage(e), call. = FALSE)
    })
    fs::file_delete(hg38_tgz)

    # download and extract mm10
    mm10_tgz <- fs::path(lib_path, "mm10_library.tar.gz")
    tryCatch({
        download.file(mm10_url, destfile = mm10_tgz, mode = "wb")
        untar(mm10_tgz, exdir = lib_path)
    }, error = function(e) {
        if (fs::file_exists(mm10_tgz)) fs::file_delete(mm10_tgz)
        stop("Failed to download or extract mm10 library: ", conditionMessage(e), call. = FALSE)
    })
    fs::file_delete(mm10_tgz)

    return(lib_path)
}

#' initialize the bart2 python package
#'
#' This function install bart2, download library and rewrite bart2 config
#'
#' This function
#'
#' @param condaenv name of conda environment that bart2 is installed
#' @param python_ver python version to use, default is "3.9"
#'
#' @return no return.
#'
#' @import magrittr reticulate fs
#'
#' @export
#'
initialize <- function(
    condaenv = NULL,
    virtualenv = NULL,
    python_ver = "3.9",
    site = "box") {
    answer1 <- askYesNo("Start installing bart2 and related data library?")

    if (is.na(answer1)) {
        message("Initiation exited")
        return(invisible(NULL))
    }

    if (!answer1) {
        message("Initiation exited")
        return(invisible(NULL))
    }

    message("Installation started...")

    # ask for existing data library path
    existing_lib <- readline("Specify the absolute path to existing data library, e.g. Data/Path/bart2_library (skip if the data library is uninstalled):")

    if (existing_lib != "") {
        if (dir.exists(existing_lib)) {
            lib_full_dir <- fs::path_norm(existing_lib)
            message("Using existing data library at: ", lib_full_dir)
        } else {
            message("Path does not exist. Proceeding to download library.")
            # download data library from selected site
            store_path <- readline("Specify the absolute path to store data library 13.3GB (skip to store under BARTsc R package directory): ")
            lib_full_dir <- get_library(store_path, site = site)
        }
    } else {
        # download data library from selected site
        store_path <- readline("Specify the absolute path to store data library 13.3GB (skip to store under BARTsc R package directory): ")
        lib_full_dir <- get_library(store_path, site = site)
    }

    message("library: ", lib_full_dir)

    install_path <- readline("Specify the absolute path to install bart2 (skip to install under BARTsc R package directory): ")
    install_bart2(install_path, python_ver = python_ver, lib = lib_full_dir)
}
