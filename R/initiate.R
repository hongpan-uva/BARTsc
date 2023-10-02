#' try to load bart2 and modules of it
#'
#' @return no return
#'
.try_import <- function() {
    tryCatch(
        {
            bart2 <<- reticulate::import("bart2", delay_load = TRUE)
            message("Load in bart2 successfully")
        },
        error = function(err) {
            stop("No available bart2", call. = FALSE)
        }
    )

    contents <- reticulate::py_get_attr(bart2, "__all__") %>% py_to_r()

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
    virtualenv = "scbart_env") {
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
#' @param lib path of bart2 library data
#'
#' @return no return
#'
#' @export
#'
#' @examples
#' install_bart2(path_of_library)
install_bart2 <- function(
    bart_dir = NULL,
    lib = NULL) {
    # create a virtual env of python 3.9:latest
    reticulate::install_python("3.9:latest")
    reticulate::virtualenv_create(envname = "scbart_env", version = "3.9:latest")
    reticulate::use_virtualenv("scbart_env")

    work_dir <- getwd()

    if (bart_dir == "") {
        bart_dir <- paste0(.libPaths()[[1]], "/scbart")
    }

    setwd(bart_dir)

    # remove existing bart2
    system("rm -rf bart2_python")

    # download latest bart2
    system2(
        command = "wget",
        args = c(
            "https://github.com/hongpan-uva/bart2/tarball/use_in_r",
            "-O", "bart2_python.tgz"
        ),
        stdout = TRUE
    )

    system("mkdir bart2_python && tar xzf bart2_python.tgz -C bart2_python --strip-components 1")
    system("rm bart2_python.tgz")
    setwd("bart2_python")

    # set up config
    # sed -i 's/=.*/= \/test1\/test2\/bart2_library/' bart.conf
    # system("sed -i 's/=.*/= \\/test1\\/test2\\/bart2_library/' bart2/bart2/bart.conf")
    lib_str <- gsub("/", "\\\\/", lib)
    paste0("sed -i 's/=.*/= ", lib_str, "/' bart2/bart.conf") %>% system()

    # install bart2 python package
    python_use <- virtualenv_python("scbart_env")
    paste(python_use, "setup.py", "install") %>% system()

    setwd(work_dir)
}

#' Download bart2 library
#'
#' @param lib_dir path to store library data
#'
#' @return path where the library date is stored
#'
#' @export
#'
#' @examples
#' get_library(path_of_library)
get_library <- function(lib_dir) {
    if (substr(lib_dir, nchar(lib_dir), nchar(lib_dir)) == "/") {
        lib_dir <- substr(lib_dir, 1, (nchar(lib_dir) - 1))
    }

    work_dir <- getwd()

    if (lib_dir == "") {
        lib_dir <- paste0(.libPaths()[[1]], "/scbart")
    }

    setwd(lib_dir)

    if ("bart2_library" %in% list.files()) {
        return(paste0(lib_dir, "/bart2_library"))
    }

    system("mkdir bart2_library")

    setwd("bart2_library")

    system2(
        command = "wget",
        args = c("https://virginia.box.com/shared/static/2kqczz9gixetcr9p4bl650uyrio5zd33.gz", "-O", "hg38_library.tar.gz"), # nolint: line_length_linter.
        stdout = TRUE
    )

    system("tar zxf hg38_library.tar.gz")

    system("rm hg38_library.tar.gz")

    system2(
        command = "wget",
        args = c("https://virginia.box.com/shared/static/bxdggnhp4bjz2l5h2zjlisnzp0ac7axf.gz", "-O", "mm10_library.tar.gz"), # nolint: line_length_linter.
        stdout = TRUE
    )

    system("tar zxf mm10_library.tar.gz")

    system("rm mm10_library.tar.gz")

    setwd(work_dir)

    return(paste0(lib_dir, "/bart2_library"))
}

#' Initiate the bart2 python package
#'
#' This function install bart2, download library and rewrite bart2 config
#'
#' This function
#'
#' @param condaenv name of conda environment that bart2 is installed
#'
#' @return no return.
#'
#' @import magrittr reticulate
#'
#' @export
#'
initiate <- function(
    condaenv = NULL,
    virtualenv = NULL) {
    answer1 <- askYesNo("Start installing bart2 and related data library?")

    if (is.na(answer1)) {
        print("Initiation exited", call. = FALSE)
        return()
    }

    if (!answer1) {
        print("Initiation exited", call. = FALSE)
        return()
    } else if (answer1) {
        message("Installation started...")

        # lib_full_dir <- "/project/zanglab_project/hz9fq/annotations/bart_library"
        lib_full_dir <- readline("Specify the path to store data library 13.3GB (skip to store under scbart R package directory): ") %>% get_library() # nolint: line_length_linter.
        print(paste("library: ", lib_full_dir))

        readline("Specify the path to install bart2 (skip to install under scbart R package directory): ") %>% install_bart2(., lib = lib_full_dir)
    }
}
