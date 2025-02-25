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
#' omnipathdb.org is ``omnipathr.license``. If you are a for-profit
#' user set it to ``'commercial'`` to make sure all the data you download
#' from OmniPath is legally allowed for commercial use. Otherwise just leave
#' it as it is: ``'academic'``.
#' If you don't use omnipathdb.org but within your organization you deployed
#' your own pypath server and want to share data whith a limited availability
#' to outside users, you may want to use a password. For this you can use
#' the ``omnipathr.password`` option.
#' Also if you want the R package to work from another pypath server instead
#' of omnipathdb.org, you can change the option ``omnipathr.url``.
#'
#' @return Nothing, this is not a function but a list.
.omnipathr_options_defaults <- list(
    omnipathr.url = 'https://omnipathdb.org/',
    omnipathr.notls_url = 'http://no-tls.omnipathdb.org/',
    omnipathr.notls_fallback = TRUE,
    omnipathr.notls_force = FALSE,
    omnipathr.license = 'academic',
    omnipathr.password = NULL,
    omnipathr.print_urls = FALSE,
    omnipathr.loglevel = 'trace',
    omnipathr.console_loglevel = 'success',
    omnipathr.logdir = NULL,
    omnipathr.logfile = NULL,
    omnipathr.cachedir = NULL,
    omnipathr.cache_timeout = 5L,
    omnipathr.use_cache = TRUE,
    omnipathr.retry_downloads = 3L,
    omnipathr.retry_downloads_in_seconds = 5L,
    omnipathr.db_lifetime = 300L,
    omnipathr.nichenet_results_dir = 'nichenet_results',
    omnipathr.uploadlists_chunk_size = 5000,
    # CURL options
    omnipathr.http_timeout = 300L,
    omnipathr.curl_verbose = FALSE,
    omnipathr.connect_timeout = 10L,
    omnipathr.tcp_keepalive = TRUE,
    omnipathr.tcp_keepintvl = 10L,
    omnipathr.tcp_keepidle = 10L,
    omnipathr.tcp_keepcnt = 13L,
    omnipathr.upkeep_interval_ms = 3e04L,
    omnipathr.ssl_verifypeer = 1L,
    omnipathr.ssl_verifyhost = 2L,
    # CURL options end
    omnipathr.uniprot_idmapping_poll_interval = 3L,
    omnipathr.uniprot_idmapping_timeout = 30L,
    omnipathr.user_agent = paste0(
        'Mozilla/5.0 (X11; Linux x86_64; rv:134.0) ',
        'Gecko/20100101 Firefox/134.0'
    ),
    omnipathr.default_http_headers = list(
        paste0(
           'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:134.0) ',
           'Gecko/20100101 Firefox/134.0'
        ),
        paste0(
           'Accept: text/html,application/xhtml+xml,',
           'application/xml;q=0.9,*/*;q=0.8'
        ),
        'Accept-Encoding: gzip, deflate, br, zstd',
        'Accept-Language: en-GB,en;q=0.5',
        'Sec-GPC: 1',
        'Upgrade-Insecure-Requests: 1',
        'Cache-Control: no-cache',
        'Pragma: no-cache'
    ),
    omnipathr.complex_translation = TRUE,
    omnipathr.complex_translation_one_to_many = 12L
)


#' Default options for a certain package
#'
#' @importFrom magrittr %>% extract2
#' @importFrom stringr str_to_lower
#' @noRd
default_options <- function(pkg = 'OmnipathR') {

    ns <- pkg %>% getNamespace

    pkg %>%
    str_to_lower %>%
    sprintf('.%s_options_defaults', .) %>%
    mget(envir = ns, ifnotfound = list(list())) %>%
    extract2(1L)

}


#' Local config file name for a certain package
#'
#' @importFrom stringr str_to_lower
#' @importFrom magrittr %>%
#' @noRd
local_config_fname <- function(pkg = 'OmnipathR') {

    pkg %>% str_to_lower %>% sprintf('%s.yml', .)

}


#' Retrieves the user level config file path of OmnipathR
#'
#' @noRd
omnipath_user_config_path <- function(){

    user_config_path('OmnipathR', 'saezlab')

}


#' Retrieves the user level config file path
#'
#' @importFrom rappdirs user_config_dir
#' @importFrom magrittr %>%
#' @importFrom stringr str_to_lower
#'
#' @noRd
user_config_path <- function(pkg = 'OmnipathR', author = 'saezlab'){

    file.path(
        user_config_dir(appname = pkg, appauthor = author),
        pkg %>% str_to_lower %>% sprintf('%s.yml', .)
    )

}


#' Default config file path of OmnipathR
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#'
#' @noRd
omnipath_default_config_path <- function(user = FALSE){

    default_config_path(user = user, pkg = 'OmnipathR', author = 'saezlab')

}


#' Default config path for a certain package
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#' @param pkg Character: name of the package
#' @param author Character: author of the package
#'
#' @return Character: path to the config file.
#'
#' @importFrom magrittr %>%
#' @noRd
default_config_path <- function(
    user = FALSE,
    pkg = 'OmnipathR',
    author = 'saezlab'
){

    pkg %>%
    local_config_fname %>%
    {`if`(file.exists(.) && !user, ., user_config_path(pkg, author))}

}


#' Current config file path of OmnipathR
#'
#' @param user Logical: prioritize the user level config even if a config in
#'     the current working directory is available.
#'
#' @return Character: path to the config file.
#'
#' @examples
#' omnipath_config_path()
#'
#' @export
omnipath_config_path <- function(user = FALSE){

    config_path(user = user, pkg = 'OmnipathR')

}


#' Current config file path for a certain package
#'
#' @param pkg Character: name of the package.
#'
#' @importFrom magrittr %>% extract2
#' @importFrom stringr str_to_lower str_to_upper
#' @export
#' @rdname omnipath_config_path
config_path <- function(user = FALSE, pkg = 'OmnipathR'){

    default <- default_config_path(user = user, pkg = pkg)
    glob_var <- pkg %>% str_to_lower %>% sprintf('%s_config', .)
    env_var <- glob_var %>% str_to_upper

    mget(
        glob_var,
        envir = .GlobalEnv,
        ifnotfound = Sys.getenv(env_var)
    ) %>%
    extract2(1L) %>%
    {`if`(nchar(.) > 0L, ., default)} %>%
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
#' options(omnipathr.license = 'commercial')
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

    save_config(path = path, title = title, local = local, pkg = 'OmnipathR')

}


#' Save the configuration of a certain package
#'
#' @param pkg Character: name of the package
#'
#' @export
#' @rdname omnipath_save_config
save_config <- function(
        path = NULL,
        title = 'default',
        local = FALSE,
        pkg = 'OmnipathR'
    ){

    path <-
        `if`(
            is.null(path),
            `if`(
                local,
                local_config_fname(pkg = pkg),
                config_path(pkg = pkg)
            ),
            path
        ) %>%
        .ensure_safe_path

    if(.path_writable(path)){

        .ensure_dir(path)
        env <- pkg %>% pkg_env

        options_to_config(pkg = pkg)
        this_config <- list()
        this_config[[title]] <- env$config

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
omnipath_load_config <- function(
        path = NULL,
        title = 'default',
        user = FALSE,
        ...
    ){

    load_config(path = path, title = title, user = user, pkg = 'OmnipathR', ...)

}


#' Load the coniguration of a certain package
#'
#' @param pkg Character: name of the package
#'
#' @importFrom yaml read_yaml write_yaml
#' @importFrom magrittr %<>%
#' @importFrom purrr map discard map_lgl
#' @export
#' @rdname omnipath_load_config
load_config <- function(
        path = NULL,
        title = 'default',
        user = FALSE,
        pkg = 'OmnipathR',
        ...
    ){

    path <-
        `if`(
            is.null(path),
            config_path(user = user, pkg = pkg),
            path
        ) %>%
        .ensure_safe_path

    yaml_config <- `if`(
        file.exists(path),
        read_yaml(file = path, ...),
        list()
    )

    if(title %in% names(yaml_config)){
        this_config <- yaml_config[[title]]
    }else{
        title <- names(yaml_config)[1L]
        if(!is.null(title)){
            this_config <- yaml_config[[title]]
        }else if(length(yaml_config) > 0L){
            this_config <- yaml_config
        }else{
            this_config <- list()
        }
    }

    if(title != 'default' && 'default' %in% names(yaml_config)){
        this_config %<>% merge_lists(yaml_config[['default']])
    }

    # formerly the parameters of this package were prefixed with "omnipath."
    # now the prefix should correspond to the package name lower case
    # here we update the existing config files
    if (this_config %>% names %>% str_detect('^omnipath\\.') %>% any) {
        names(this_config) %<>% str_replace('^omnipath\\.', 'omnipathr.')
        write_yaml(list(this_config) %>% set_names(title), file = path)
    }

    pkg_l <- pkg %>% pkg_prefix('', pkg = .)

    # ensure the `pkg_lower.` prefix for all parameter keys
    names(this_config) <- ifelse(
        startsWith(as.character(names(this_config)), pkg_l),
        names(this_config),
        sprintf('%s%s', pkg_l, names(this_config))
    )

    this_config %<>% merge_lists(default_options(pkg))
    this_config %<>% merge_lists(env_config(pkg), .)

    env <- pkg %>% pkg_env
    env$config <- this_config

    options_to_config(pkg = pkg)
    config_to_options(pkg = pkg)

}


#' Options from environment variables
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom purrr map2 discard
#' @importFrom stringr str_replace_all str_to_upper
#' @noRd
env_config <- function(pkg = 'OmnipathR'){

    pkg %>%
    default_options %>%
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


#' Restore the built-in default values of all config parameters of a package
#'
#' @param pkg Character: name of a package
#'
#' @importFrom magrittr %>%
#' @export
#' @rdname omnipath_reset_config
reset_config <- function(save = NULL, reset_all = FALSE, pkg = 'OmnipathR'){

    env <- pkg %>% pkg_env
    env$config <- pkg %>% default_options

    if(!reset_all){
        options_to_config(pkg = pkg)
    }

    config_to_options(pkg = pkg)

    if(!is.null(save)){
        path <- `if`(
            is.logical(save) && save,
            config_path(pkg = pkg),
            save
        )
        save_config(path, pkg = pkg)
    }

    invisible(env$config)

}


#' Restore the built-in default values of all config parameters of OmnipathR
#'
#' @param save If a path, the restored config will be also saved
#'     to this file. If TRUE, the config will be saved to the current default
#'     config path (see \code{\link{omnipath_config_path}}).
#' @param reset_all Reset to their defaults also the options already set in
#'     the R options.
#' @param ... Ignored.
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
#' @importFrom purrr partial
#' @seealso \code{\link{omnipath_load_config}, \link{omnipath_save_config}}
omnipath_reset_config <- partial(reset_config, pkg = 'OmnipathR')


#' Populates the config from the default local or user level config file
#' or the built-in defaults.
#'
#' @noRd
init_config <- function(user = FALSE, pkg = 'OmnipathR'){

    config_path(user = user, pkg = pkg) %>%
    {`if`(
        file.exists(.),
        load_config(., pkg = pkg),
        # this normally happens only at the very first use of the package
        reset_config(save = ., pkg = pkg)
    )}

}


#' Populate the config of OmnipathR
#'
#' @importFrom purrr partial
#' @noRd
omnipath_init_config <- partial(init_config, pkg = 'OmnipathR')


#' Loads the settings from omnipathr.env$config to options
#'
#' @noRd
config_to_options <- function(pkg = 'OmnipathR'){

    env <- pkg %>% pkg_env
    do.call(options, env$config)

}


#' Copy a package's settings from options to env$config
#'
#' @importFrom magrittr %>%
#' @importFrom purrr discard
#' @noRd
options_to_config <- function(pkg = 'OmnipathR'){

    env <- pkg %>% pkg_env

    env$config <-
        from_options <- do.call(
            options,
            as.list(names(default_options(pkg = pkg)))
        ) %>%
        discard(is.null) %>%
        merge_lists(env$config)

}


#' Populates the URL register
#'
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @noRd
.load_urls <- function(pkgname){

    env <- pkgname %>% pkg_env

    env$urls <-
        system.file(
            'internal',
            'urls.json',
            package = pkgname,
            mustWork = TRUE
        ) %>%
        safe_json(simplifyVector = TRUE) %>%
        map(paste0, collapse = '')

}


#' Set up the logfile and logging parameters for OmnipathR
#'
#' @noRd
omnipath_init_log <- function(){

    init_log('OmnipathR')

}


#' Set up the logfile and logging parameters.
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom rlang %||%
#' @importFrom logger log_formatter log_appender log_layout appender_file
#' @importFrom logger appender_console log_threshold layout_glue_generator
#' @importFrom logger formatter_glue_or_sprintf
#' @noRd
init_log <- function(pkg = 'OmnipathR'){

    logfile_enabled <- pkg %>% has_logfile()
    pkg_l <- pkg %>% str_to_lower

    if(logfile_enabled){

        log_path <-
            'logfile' %>%
            pkg_getopt(pkg) %||%
            file.path(
                pkg_getopt('logdir', pkg) %||% sprintf('%s-log', pkg_l),
                sprintf('%s-%s.log', pkg_l, format(Sys.time(), "%Y%m%d-%H%M"))
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
            '%sloglevel',
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
            namespace = pkg,
            index = idx
        )
        log_threshold(
            ensure_loglevel(pkg_getopt(loglevel, pkg)),
            namespace = pkg,
            index = idx
        )
        log_appender(appender, namespace = pkg, index = idx)
        log_layout(
            layout_glue_generator(format = layout_format),
            namespace = pkg,
            index = idx
        )

    }

}
