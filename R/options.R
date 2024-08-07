#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2024
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Alberto Valdeolivas
#                  Dénes Türei (turei.denes@gmail.com)
#                  Attila Gábor
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Website: https://r.omnipathdb.org/
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Default values for the package options
#'
#' These options describe the default settings for OmnipathR so you do not
#' need to pass these parameters at each function call.
#' Currently the only option useful for the public web service at
#' omnipathdb.org is ``omnipath.license``. If you are a for-profit
#' user set it to ``'commercial'`` to make sure all the data you download
#' from OmniPath is legally allowed for commercial use. Otherwise just leave
#' it as it is: ``'academic'``.
#' If you don't use omnipathdb.org but within your organization you deployed
#' your own pypath server and want to share data whith a limited availability
#' to outside users, you may want to use a password. For this you can use
#' the ``omnipath.password`` option.
#' Also if you want the R package to work from another pypath server instead
#' of omnipathdb.org, you can change the option ``omnipath.url``.
#'
#' @return Nothing, this is not a function but a list.
.omnipath_options_defaults <- list(
    omnipath.url = 'https://omnipathdb.org/',
    omnipath.notls_url = 'http://no-tls.omnipathdb.org/',
    omnipath.notls_fallback = TRUE,
    omnipath.notls_force = FALSE,
    omnipath.license = 'academic',
    omnipath.password = NULL,
    omnipath.print_urls = FALSE,
    omnipath.loglevel = 'trace',
    omnipath.console_loglevel = 'success',
    omnipath.logdir = NULL,
    omnipath.logfile = NULL,
    omnipath.cachedir = NULL,
    omnipath.cache_timeout = 5,
    omnipath.use_cache = TRUE,
    omnipath.retry_downloads = 3,
    omnipath.retry_downloads_in_seconds = 5,
    omnipath.db_lifetime = 300,
    omnipath.nichenet_results_dir = 'nichenet_results',
    omnipath.uploadlists_chunk_size = 5000,
    omnipath.connect_timeout = 10,
    omnipath.uniprot_idmapping_poll_interval = 3,
    omnipath.uniprot_idmapping_timeout = 30,
    omnipath.user_agent = paste0(
        'Mozilla/5.0 (X11; Linux x86_64; rv:120.0) ',
        'Gecko/20100101 Firefox/120.0'
    ),
    omnipath.complex_translation = TRUE,
    omnipath.complex_translation_one_to_many = TRUE
)


.omnipath_local_config_fname <- 'omnipathr.yml'


#' Retrieves the user level config file path
#'
#' @importFrom rappdirs user_config_dir
#'
#' @noRd
omnipath_get_user_config_path <- function(){

    file.path(
        user_config_dir(appname = 'OmnipathR', appauthor = 'saezlab'),
        'omnipathr.yml'
    )

}


#' Retrieves the default config file path
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#'
#' @noRd
omnipath_get_default_config_path <- function(user = FALSE){

    `if`(
        file.exists(.omnipath_local_config_fname) && !user,
        .omnipath_local_config_fname,
        omnipath_get_user_config_path()
    )

}


#' Current config file path
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#'
#' @return Character: path to the config file.
#'
#' @examples
#' omnipath_get_config_path()
#'
#' @importFrom magrittr %>%
#' @export
omnipath_get_config_path <- function(user = FALSE){

    config_path_default <- omnipath_get_default_config_path(user = user)

    config_path <- mget(
        'omnipath_config',
        envir = .GlobalEnv,
        ifnotfound = Sys.getenv('OMNIPATHR_CONFIG')
    )[[1]]

    `if`(
        nchar(config_path) > 0,
        config_path,
        config_path_default
    ) %>%
    .ensure_safe_path

}


#' Save the current package configuration
#'
#' @param path Path to the config file. Directories and the file will be
#'     created if don't exist.
#' @param title Save the config under this title. One config file might
#'     contain multiple configurations, each identified by a title.
#' @param local Save into a config file in the current directory instead of
#'     a user level config file. When loading, the config in the current
#'     directory has priority over the user level config.
#'
#' @return Returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # after this, all downloads will default to commercial licenses
#' # i.e. the resources that allow only academic use will be excluded:
#' options(omnipath.license = 'commercial')
#' omnipath_save_config()
#' }
#'
#' @export
#' @importFrom yaml write_yaml
omnipath_save_config <- function(
        path = NULL,
        title = 'default',
        local = FALSE
    ){

    path <-
        `if`(
            is.null(path),
            `if`(
                local,
                .omnipath_local_config_fname,
                omnipath_get_config_path()
            ),
            path
        ) %>%
        .ensure_safe_path

    if(.path_writable(path)){

        .ensure_dir(path)

        omnipath_options_to_config()
        this_config <- list()
        this_config[[title]] <- omnipath.env$config

        write_yaml(this_config, file = path)

    }

}


#' Load the package configuration from a config file
#'
#' @param path Path to the config file.
#' @param title Load the config under this title. One config file might
#'     contain multple configurations, each identified by a title. If the
#'     title is not available the first section of the config file will be
#'     used.
#' @param user Force to use the user level config even if a config file
#'     exists in the current directory. By default, the local config files
#'     have prioroty over the user level config.
#' @param ... Passed to \code{yaml::yaml.load_file}.
#'
#' @return Invisibly returns the config as a list.
#'
#' @examples
#' \dontrun{
#' # load the config from a custom config file:
#' omnipath_load_config(path = 'my_custom_omnipath_config.yml')
#' }
#'
#' @export
#' @importFrom yaml read_yaml write_yaml
#' @importFrom magrittr %<>%
#' @importFrom purrr map discard map_lgl
omnipath_load_config <- function(
        path = NULL,
        title = 'default',
        user = FALSE,
        ...
    ){

    path <-
        `if`(
            is.null(path),
            omnipath_get_config_path(user = user),
            path
        ) %>%
        .ensure_safe_path

    yaml_config <- `if`(
        file.exists(path),
        read_yaml(file = path, ...),
        list()
    )

    # Earlier we had the URLs in the config. Then I realized it's not a
    # good practice, and apart from not having them any more, here we
    # remove them from existing configs.
    if(.path_writable(path) && yaml_config %>% map_lgl(is.list) %>% all){

        yaml_config %<>% map(
            ~`[`(
                .x,
                .x %>% names %>% discard(endsWith, '_url')
            )
        )
        write_yaml(yaml_config, file = path)

    }

    if(title %in% names(yaml_config)){
        this_config <- yaml_config[[title]]
    }else{
        title <- names(yaml_config)[1]
        if(!is.null(title)){
            this_config <- yaml_config[[title]]
        }else if(length(yaml_config) > 0){
            this_config <- yaml_config
        }else{
            this_config <- list()
        }
    }

    if(title != 'default' && 'default' %in% names(yaml_config)){
        this_config %<>% merge_lists(yaml_config[['default']])
    }

    # ensure the `omnipath.` prefix for all parameter keys
    names(this_config) <- ifelse(
        startsWith(as.character(names(this_config)), 'omnipath.'),
        names(this_config),
        sprintf('omnipath.%s', names(this_config))
    )

    this_config %<>% merge_lists(.omnipath_options_defaults)
    this_config %<>% merge_lists(omnipath_env_config(), .)

    omnipath.env$config <- this_config

    omnipath_options_to_config()
    omnipath_config_to_options()

}


#' Options from environment variables
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map2 discard
#' @importFrom stringr str_replace_all str_to_upper
#' @noRd
omnipath_env_config <- function(){

    .omnipath_options_defaults %>%
    map2(
        names(.),
        function(val, key){
            key %>%
            str_replace_all('\\.', '_') %>%
            str_to_upper %>%
            Sys.getenv %>%
            as_type(val)
        }
    ) %>%
    discard(function(x){is.na(x) || x == ''})

}


#' Restores the built-in default values of all config parameters
#'
#' @param save If a path, the restored config will be also saved
#'     to this file. If TRUE, the config will be saved to the current default
#'     config path (see \code{\link{omnipath_get_config_path}}).
#' @param reset_all Reset to their defaults also the options already set in
#'     the R options.
#'
#' @examples
#' \dontrun{
#' # restore the defaults and write them to the default config file:
#' omnipath_reset_config()
#' omnipath_save_config()
#' }
#'
#' @return The config as a list.
#'
#' @export
#' @seealso \code{\link{omnipath_load_config}, \link{omnipath_save_config}}
omnipath_reset_config <- function(save = NULL, reset_all = FALSE){

    omnipath.env$config <- .omnipath_options_defaults

    if(!reset_all){
        omnipath_options_to_config()
    }

    omnipath_config_to_options()

    if(!is.null(save)){
        path <- `if`(
            is.logical(save) && save,
            omnipath_get_config_path(),
            save
        )
        omnipath_save_config(path)
    }

    invisible(omnipath.env$config)

}


#' Populates the config from the default local or user level config file
#' or the built-in defaults.
#'
#' @noRd
omnipath_init_config <- function(user = FALSE){

    config_path <- omnipath_get_config_path(user = user)

    if(file.exists(config_path)){
        omnipath_load_config(config_path)
    }else{
        # this normally happens only at the very first use of the package
        omnipath_reset_config(save = config_path)
    }

}


#' Loads the settings from omnipath.env$config to options
#'
#' @noRd
omnipath_config_to_options <- function(){

    do.call(options, omnipath.env$config)

}

#' Copies OmnipathR settings from options to omnipath.env$config
#'
#' @importFrom magrittr %>%
#' @importFrom purrr discard
#' @noRd
omnipath_options_to_config <- function(){

    from_options <- do.call(
        options,
        as.list(names(.omnipath_options_defaults))
    ) %>%
    discard(is.null)

    omnipath.env$config <- merge_lists(
        from_options,
        omnipath.env$config
    )

}


#' Populates the URL register
#'
#' @importFrom purrr map
#'
#' @noRd
.load_urls <- function(pkgname){

    omnipath.env$urls <-
        system.file(
            'internal',
            'urls.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        safe_json(simplifyVector = TRUE) %>%
        map(paste0, collapse = '')

}


#' Setting up the logfile and logging parameters.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang %||%
#' @importFrom logger log_formatter log_appender log_layout appender_file
#' @importFrom logger appender_console log_threshold layout_glue_generator
#' @importFrom logger formatter_glue_or_sprintf
#'
#' @noRd
omnipath_init_log <- function(pkgname = 'OmnipathR'){

    logfile_enabled <- omnipath_has_logfile()

    if(logfile_enabled){

        log_path <-
            getOption('omnipath.logfile') %||%
            file.path(
                getOption('omnipath.logdir') %||% 'omnipathr-log',
                sprintf('omnipathr-%s.log', format(Sys.time(), "%Y%m%d-%H%M"))
            ) %>%
            absolute_path() %T>%
            {dir.create(
                dirname(.),
                showWarnings = FALSE,
                recursive = TRUE
            )} %>%
            .ensure_safe_path

    }

    for(idx in seq(1 + logfile_enabled)){
        # 1 = console
        # 2 = logfile

        loglevel <- sprintf(
            'omnipath.%sloglevel',
            `if`(idx == 1, 'console_', '')
        )
        appender <- `if`(
            idx == 1,
            logger::appender_console,
            logger::appender_file(log_path)
        )
        layout_format <- `if`(
            idx == 1,
            paste0(
                '[{format(time, "%Y-%m-%d %H:%M:%S")}] ',
                '[{colorize_by_log_level(level, levelr)}]',
                '{paste0(rep(" ", 7 - nchar(level)), collapse = "")} ',
                '[{ns}] ',
                '{grayscale_by_log_level(msg, levelr)}'
            ),
            paste0(
                '[{format(time, "%Y-%m-%d %H:%M:%S")}] ',
                '[{level}]',
                '{paste0(rep(" ", 7 - nchar(level)), collapse = "")} ',
                '[{ns}] ',
                '{msg}'
            )
        )

        log_formatter(
            formatter_glue_or_sprintf,
            namespace = pkgname,
            index = idx
        )
        log_threshold(
            ensure_loglevel(options(loglevel)[[1]]),
            namespace = pkgname,
            index = idx
        )
        log_appender(appender, namespace = pkgname, index = idx)
        log_layout(
            layout_glue_generator(format = layout_format),
            namespace = pkgname,
            index = idx
        )

    }

}
