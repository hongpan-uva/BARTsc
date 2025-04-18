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
    out <- tryCatch(
        {
            reticulate::use_python_version(python_ver)
        },
        error = function(cond) {
            message("Given python version was not found, will install it")
            return(1)
        }
    )

    # install certain python version
    if (out == 1) {
        reticulate::install_python(python_ver)
    }

    # create a virtual env
    reticulate::virtualenv_create(envname = "bartsc_env", version = python_ver)
    reticulate::use_virtualenv("bartsc_env")
    reticulate::virtualenv_remove(envname = "bartsc_env", packages = "numpy", confirm = FALSE) # remove existing numpy

    work_dir <- getwd()

    if (is.null(bart_dir) || bart_dir == "") {
        bart_dir <- paste0(.libPaths()[[1]], "/BARTsc")
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
    python_use <- reticulate::virtualenv_python("bartsc_env")
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
get_library <- function(lib_dir) {
    if (substr(lib_dir, nchar(lib_dir), nchar(lib_dir)) == "/") {
        lib_dir <- substr(lib_dir, 1, (nchar(lib_dir) - 1))
    }

    work_dir <- getwd()

    if (lib_dir == "") {
        lib_dir <- paste0(.libPaths()[[1]], "/BARTsc")
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
#' @import magrittr reticulate
#'
#' @export
#'
initialize <- function(
    condaenv = NULL,
    virtualenv = NULL,
    python_ver = "3.9") {
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
        lib_full_dir <- readline("Specify the path to store data library 13.3GB (skip to store under BARTsc R package directory): ") %>% get_library() # nolint: line_length_linter.
        print(paste("library: ", lib_full_dir))

        readline("Specify the path to install bart2 (skip to install under BARTsc R package directory): ") %>% install_bart2(., python_ver = python_ver, lib = lib_full_dir)
    }
}
