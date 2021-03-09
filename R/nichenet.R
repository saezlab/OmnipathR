#!/usr/bin/env Rscript

#
#  This file is part of the `OmnipathR` R package
#
#  Copyright
#  2018-2021
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
#  Website: https://saezlab.github.io/omnipathr
#  Git repo: https://github.com/saezlab/OmnipathR
#


#' Executes the full NicheNet pipeline
#'
#' Builds all prior knowledge data required by NicheNet. For this it calls
#' a multitude of methods to download and combine data from various
#' databases according to the settings. The content of the prior knowledge
#' data is highly customizable, see the documentation of the related
#' functions. After the prior knowledge is ready, it performs parameter
#' optimization to build a NicheNet model. This results a weighted ligand-
#' target matrix. Then, considering the expressed genes from user provided
#' data, a gene set of interest and background genes, it executes the
#' NicheNet ligand activity analysis.
#'
#' @param only_omnipath Logical: use only OmniPath for network knowledge.
#'     This is a simple switch for convenience, further options are available
#'     by the other arguments. By default we use all available resources.
#'     The networks can be customized on a resource by resource basis, as
#'     well as providing custom parameters for individual resources, using
#'     the parameters `signaling_network`, `lr_network` and `gr_network`.
#' @param expressed_genes_transmitter Character vector with the gene symbols
#'     of the genes expressed in the cells transmitting the signal.
#' @param expressed_genes_receiver Character vector with the gene symbols
#'     of the genes expressed in the cells receiving the signal.
#' @param genes_of_interest Character vector with the gene symbols of the
#'     genes of interest. These are the genes in the receiver cell population
#'     that are potentially affected by ligands expressed by interacting
#'     cells (e.g. genes differentially expressed upon cell-cell interaction).
#' @param background_genes Character vector with the gene symbols of the
#'     genes to be used as background.
#' @param n_top_ligands How many of the top ligands to include in the
#'     ligand-target table.
#' @param n_top_targets How many of the top targets (for each of the top
#'     ligands) to consider in the ligand-target table.
#' @param signaling_network A list of parameters for building the signaling
#'     network, passed to \code{\link{nichenet_signaling_network}}
#' @param lr_network A list of parameters for building the ligand-receptor
#'     network, passed to \code{\link{nichenet_lr_network}}
#' @param gr_network A list of parameters for building the gene regulatory
#'     network, passed to \code{\link{nichenet_gr_network}}
#' @param make_multi_objective_function_param Override parameters for
#'     \code{smoof::makeMultiObjectiveFunction}.
#' @param objective_function_param Override additional arguments passed to
#'     the objective function.
#' @param mlrmbo_optimization_param Override arguments for
#'     \code{nichenetr::mlrmbo_optimization}.
#' @param construct_ligand_target_matrix_param Override parameters for
#'     \code{nichenetr::construct_ligand_target_matrix}.
#' @param results_dir Character: path to the directory to save intermediate
#'     and final outputs from NicheNet methods.
#'
#' @examples
#' \donttest{
#' nichenet_results <- nichenet_main(
#'     # altering some network resource parameters, the rest
#'     # of the resources will be loaded according to the defaults
#'     signaling_network = list(
#'         cpdb = NULL, # this resource will be excluded
#'         evex = list(min_confidence = 1.0) # override some parameters
#'     ),
#'     gr_network = list(only_omnipath = TRUE),
#'     n_top_ligands = 20,
#'     # override the default number of CPU cores to use
#'     mlrbo_optimization_param = list(ncores = 4)
#' )
#' }
#'
#' @export
#' @importFrom magrittr %>% %T>%
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_networks}}}
#'     \item{\code{\link{nichenet_signaling_network}}}
#'     \item{\code{\link{nichenet_lr_network}}}
#'     \item{\code{\link{nichenet_gr_network}}}
#' }
nichenet_main <- function(
    only_omnipath = FALSE,
    expressed_genes_transmitter = NULL,
    expressed_genes_receiver = NULL,
    genes_of_interest = NULL,
    background_genes = NULL,
    n_top_ligands = 42,
    n_top_targets = 250,
    signaling_network = list(),
    lr_network = list(),
    gr_network = list(),
    make_multi_objective_function_param = list(),
    objective_function_param = list(),
    mlrmbo_optimization_param = list(),
    construct_ligand_target_matrix_param = list(),
    results_dir = NULL
){

    # NSE vs. R CMD check workaround
    optimization_results <- optimized_parameters <- weighted_networks <-
        ligand_target_matrix <- NULL

    top_env <- environment()

    results_dir %>%
    if_null(options('omnipath.nichenet_results_dir')) %>%
    options(omnipath.nichenet_results_dir = .)

    logger::log_success('Building NicheNet prior knowledge.')

    networks <- nichenet_networks(
        signaling_network = signaling_network,
        lr_network = lr_network,
        gr_network = gr_network,
        only_omnipath = only_omnipath
    )

    expression <- nichenet_expression_data() %>%
        nichenet_remove_orphan_ligands(lr_network = networks$lr_network)

    logger::log_success('Finished building NicheNet prior knowledge.')
    logger::log_success('Building NicheNet model.')

    networks %>%
    nichenet_optimization(
        expression = expression,
        make_multi_objective_function_param =
            make_multi_objective_function_param,
        objective_function_param = objective_function_param,
        mlrmbo_optimization_param = mlrmbo_optimization_param
    ) %T>%
    assign(x = 'optimization_results', value = ., envir = top_env) %>%
    nichenet_build_model(networks = networks, weighted = FALSE) %>%
    c(
        list(
            lr_network = networks$lr_network,
            construct_ligand_target_matrix_param =
                construct_ligand_target_matrix_param,
            weighted = FALSE
        )
    ) %T>%
    assign(
        x = 'optimized_parameters',
        value = .$optimized_parameters,
        envir = top_env
    ) %>%
    do.call(what = nichenet_ligand_target_matrix) %>%
    # second run with weights
    {nichenet_build_model(
        networks = networks,
        optimization_results = optimization_results
    )} %T>%
    assign(
        x = 'weighted_networks',
        value = .$weighted_networks,
        envir = top_env
    ) %>%
    c(
        list(
            lr_network = networks$lr_network,
            construct_ligand_target_matrix_param =
                construct_ligand_target_matrix_param
        )
    ) %>%
    do.call(what = nichenet_ligand_target_matrix) %T>%
    assign(x = 'ligand_target_matrix', ., envir = top_env) %T>%
    {logger::log_success('Finished building NicheNet model.')} %>%
    list(
        ligand_target_matrix = .,
        expressed_genes_transmitter = expressed_genes_transmitter,
        expressed_genes_receiver = expressed_genes_receiver,
        genes_of_interest = genes_of_interest,
        background_genes = background_genes
    ) %>%
    {`if`(
        any(map_lgl(., is.null)),
        list(ligand_activities = NULL, ligand_target_links = NULL),
        do.call(nichenet_ligand_activities, .)
    )} %>%
    c(
        list(
            networks = networks,
            expression = expression,
            optimized_parameters = optimized_parameters,
            weighted_networks = weighted_networks,
            ligand_target_matrix = ligand_target_matrix
        ),
        .
    ) %T>%
    {logger::log_success('Completed NicheNet pipeline.')}

}


#' Removes experiments with orphan ligands
#'
#' Removes from the expression data the perturbation experiments involving
#' ligands without connections.
#'
#' @param expression Expression data as returned by
#'     \code{\link{nichenet_expression_data}}.
#' @param lr_network A NicheNet format ligand-recptor network data frame as
#'     produced by \code{\link{nichenet_lr_network}}.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' expression <- nichenet_remove_orphan_ligands(
#'     expression,
#'     networks$lr_network
#' )
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom purrr keep
nichenet_remove_orphan_ligands <- function(expression, lr_network){

    # NSE vs. R CMD check workaround
    networks <- NULL

    all_ligands <- networks$lr_network$from %>% unique

    expression %>%
    keep(
        function(record){
            all(record$from %in% all_ligands)
        }
    )

}


#' Optimizes NicheNet model parameters
#'
#' Optimize NicheNet method parameters, i.e. PageRank parameters and source
#' weights, basedon a collection of experiments where the effect of a ligand
#' on gene expression was measured.
#'
#' @param networks A list with NicheNet format signaling, ligand-receptor
#'     and gene regulatory networks as produced by
#'     \code{\link{nichenet_networks}}.
#' @param expression A list with expression data from ligand perturbation
#'     experiments, as produced by \code{\link{nichenet_expression_data}}.
#' @param make_multi_objective_function_param Override parameters for
#'     \code{smoof::makeMultiObjectiveFunction}.
#' @param objective_function_param Override additional arguments passed to
#'     the objective function.
#' @param mlrmbo_optimization_param Override arguments for
#'     \code{nichenetr::mlrmbo_optimization}.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' optimization_results <- nichenet_optimization(networks, expression)
#' }
#'
#' @export
#' @importFrom magrittr %>% %T>% %<>%
#' @importFrom purrr map
#' @importFrom dplyr mutate
nichenet_optimization <- function(
    networks,
    expression,
    make_multi_objective_function_param = list(),
    objective_function_param = list(),
    mlrmbo_optimization_param = list()
){

    resources <-
        networks %>%
        map(list('source', unique)) %>%
        unlist %>%
        unique

    n_resources <- resources %>% length

    optimization_results_rds_path <-
        nichenet_results_dir() %>%
        file.path('optimization_results.rds')

    mof_topology_correction <-
        make_multi_objective_function_param %>%
        add_defaults(
            fun = smoof::makeMultiObjectiveFunction,
            defaults = list(
                name = 'nichenet_optimization',
                description = paste(
                    'Data source weight and hyperparameter optimization:',
                    'expensive black-box function'
                ),
                fn = nichenetr::model_evaluation_optimization,
                par.set = ParamHelpers::makeParamSet(
                    ParamHelpers::makeNumericVectorParam(
                        'source_weights',
                        len = n_resources,
                        lower = 0,
                        upper = 1,
                        tunable = FALSE
                    ),
                    ParamHelpers::makeNumericVectorParam(
                        'lr_sig_hub',
                        len = 1,
                        lower = 0,
                        upper = 1,
                        tunable = TRUE
                    ),
                    ParamHelpers::makeNumericVectorParam(
                        'gr_hub',
                        len = 1,
                        lower = 0,
                        upper = 1,
                        tunable = TRUE
                    ),
                    ParamHelpers::makeNumericVectorParam(
                        'ltf_cutoff',
                        len = 1,
                        lower = 0.9,
                        upper = 0.999,
                        tunable = TRUE
                    ),
                    ParamHelpers::makeNumericVectorParam(
                        'damping_factor',
                        len = 1,
                        lower = 0.01,
                        upper = 0.99,
                        tunable =TRUE
                    )
                ),
                has.simple.signature = FALSE,
                n.objectives = 4,
                noisy = FALSE,
                minimize = rep(FALSE, 4)
            )
        ) %>%
        do.call(what = smoof::makeMultiObjectiveFunction)

    objective_function_param %<>%
        add_defaults(
            fun = mof_topology_correction,
            defaults = list(
                source_names = resources,
                algorithm = 'PPR',
                correct_topology = FALSE,
                lr_network = networks$lr_network,
                sig_network = networks$signaling_network,
                gr_network = networks$gr_network,
                settings = expression %>% map(
                    nichenetr::convert_expression_settings_evaluation
                ),
                secondary_targets = FALSE,
                remove_direct_links = 'no',
                cutoff_method = 'quantile'
            )
        )

    mlrmbo_optimization_param %>%
    add_defaults(
        fun = nichenetr::mlrmbo_optimization,
        defaults = list(
            run_id = 1,
            obj_fun = mof_topology_correction,
            niter = 8,
            ncores = 8,
            nstart = 160,
            additional_arguments = objective_function_param
        )
    ) %T>%
    {logger::log_success(
        'Running multi-objective model-based optimization'
    )} %>%
    do.call(what = nichenetr::mlrmbo_optimization) %T>%
    saveRDS(optimization_results_rds_path) %T>%
    {logger::log_success(
        paste0(
            'Multi-objective model-based optimization ready, ',
            'saved results to `%s`.'
        ),
        optimization_results_rds_path
    )}

}


#' Construct a NicheNet ligand-target model
#'
#' @param optimization_results The outcome of NicheNet parameter optimization
#'     as produced by \code{\link{nichenet_optimization}}.
#' @param networks A list with NicheNet format signaling, ligand-receptor
#'     and gene regulatory networks as produced by
#'     \code{\link{nichenet_networks}}.
#' @param weighted Logical: whether to use the optimized weights.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' optimization_results <- nichenet_optimization(networks, expression)
#' nichenet_model <- nichenet_build_model(optimization_results, networks)
#' }
#'
#' @export
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom magrittr %>% %T>%
#' @importFrom dplyr pull
nichenet_build_model <- function(
    optimization_results,
    networks,
    weighted = TRUE
){

    # all resources with initial weights of 1
    resource_weights <-
        networks %>%
        map(list('source', unique)) %>%
        unlist %>%
        unique %>%
        {tibble(source = .)} %>%
        mutate(weight = 1)

    optimized_parameters <-
        optimization_results %T>%
        {logger::log_success('Processing MLRMBO parameters.')} %>%
        nichenetr::process_mlrmbo_nichenet_optimization(
            source_names = resource_weights %>% pull(source) %>% unique
        ) %T>%
        {logger::log_success('Finished processing MLRMBO parameters.')}

    weighted_networks_rds_path <-
        nichenet_results_dir() %>%
        file.path(
            sprintf(
                'weighted_networks%s.rds',
                `if`(weighted, '_weighted', '')
            )
        )

    logger::log_success('Creating weighted networks.')

    nichenetr::construct_weighted_networks(
        lr_network = networks$lr_network,
        sig_network = networks$signaling_network,
        gr_network = networks$gr_network,
        source_weights_df = `if`(
            weighted,
            optimization_results$source_weight_df,
            resource_weights
        )
    ) %T>%
    {logger::log_success('Applying hub corrections.')} %>%
    nichenetr::apply_hub_corrections(
        lr_sig_hub = optimized_parameters$lr_sig_hub,
        gr_hub = optimized_parameters$gr_hub
    ) %T>%
    saveRDS(weighted_networks_rds_path) %T>%
    {logger::log_success(
        'Created weighted networks, saved to `%s`.',
        weighted_networks_rds_path
    )} %>%
    list(
        weighted_networks = .,
        optimized_parameters = optimized_parameters
    )

}


#' Creates a NicheNet ligand-target matrix
#'
#' @param weighted_networks Weighted networks as pro
#' @param lr_network A data frame with ligand-receptor interactions, as
#'     produced by \code{\link{nichenet_lr_network}}.vided by
#'     \code{\link{nichenet_build_model}}
#' @param optimized_parameters The outcome of NicheNet parameter optimization
#'     as produced by \code{\link{nichenet_build_model}}.
#' @param weighted Logical: wether the network sources are weighted. In this
#'     function it only affects the output file name.
#' @param construct_ligand_target_matrix_param Override parameters for
#'     \code{nichenetr::construct_ligand_target_matrix}.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' optimization_results <- nichenet_optimization(networks, expression)
#' nichenet_model <- nichenet_build_model(optimization_results, networks)
#' lt_matrix <- nichenet_ligand_target_matrix(
#'     nichenet_model$weighted_networks,
#'     networks$lr_network,
#'     nichenet_model$optimized_parameters
#' )
#' }
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom magrittr %>% %T>%
nichenet_ligand_target_matrix <- function(
    weighted_networks,
    lr_network,
    optimized_parameters,
    weighted = TRUE,
    construct_ligand_target_matrix_param = list()
){

    # NSE vs. R CMD check workaround
    from <- NULL

    ligands <- lr_network %>% pull(from) %>% unique %>% as.list

    ligand_target_matrix_rds_path <-
        nichenet_results_dir() %>%
        file.path(
            sprintf(
                'ligand_target_matrix%s.rds',
                `if`(weighted, '_weighted', '')
            )
        )

    construct_ligand_target_matrix_param %>%
    add_defaults(
        fun = nichenetr::construct_ligand_target_matrix,
        defaults = list(
            weighted_networks = weighted_networks,
            ligands = ligands,
            algorithm = 'PPR',
            damping_factor = optimized_parameters$damping_factor,
            ltf_cutoff = optimized_parameters$ltf_cutoff
        )
    ) %T>%
    {logger::log_success('Creating ligand-target matrix.')} %>%
    do.call(what = nichenetr::construct_ligand_target_matrix) %T>%
    saveRDS(ligand_target_matrix_rds_path) %T>%
    {logger::log_success(
        'Created ligand-target matrix, saved to `%s`.',
        ligand_target_matrix_rds_path
    )}

}


#' Calls the NicheNet ligand activity analysis
#'
#' @param ligand_target_matrix A matrix with rows and columns corresponding
#'     to ligands and targets, respectively. Produced by
#'     \code{\link{nichenet_ligand_target_matrix}} or
#'     \code{nichenetr::construct_ligand_target_matrix}.
#' @param lr_network A data frame with ligand-receptor interactions, as
#'     produced by \code{\link{nichenet_lr_network}}.
#' @param expressed_genes_transmitter Character vector with the gene symbols
#'     of the genes expressed in the cells transmitting the signal.
#' @param expressed_genes_receiver Character vector with the gene symbols
#'     of the genes expressed in the cells receiving the signal.
#' @param genes_of_interest Character vector with the gene symbols of the
#'     genes of interest. These are the genes in the receiver cell population
#'     that are potentially affected by ligands expressed by interacting
#'     cells (e.g. genes differentially expressed upon cell-cell interaction).
#' @param background_genes Character vector with the gene symbols of the
#'     genes to be used as background.
#' @param n_top_ligands How many of the top ligands to include in the
#'     ligand-target table.
#' @param n_top_targets For each ligand, how many of the top targets to
#'     include in the ligand-target table.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' optimization_results <- nichenet_optimization(networks, expression)
#' nichenet_model <- nichenet_build_model(optimization_results, networks)
#' lt_matrix <- nichenet_ligand_target_matrix(
#'     nichenet_model$weighted_networks,
#'     networks$lr_network,
#'     nichenet_model$optimized_parameters
#' )
#' ligand_activities <- nichenet_ligand_activities(
#'     ligand_target_matrix = lt_matrix,
#'     lr_network = networks$lr_network,
#'     # the rest of the parameters should come
#'     # from your transcriptomics data:
#'     expressed_genes_transmitter = expressed_genes_transmitter,
#'     expressed_genes_receiver = expressed_genes_receiver,
#'     genes_of_interest = genes_of_interest
#' )
#' }
#'
#' @export
#' @importFrom magrittr %>% %<>% %T>%
#' @importFrom dplyr pull filter arrange
nichenet_ligand_activities <- function(
    ligand_target_matrix,
    lr_network,
    expressed_genes_transmitter,
    expressed_genes_receiver,
    genes_of_interest,
    background_genes = NULL,
    n_top_ligands = 42,
    n_top_targets = 250
){

    # NSE vs. R CMD check workaround
    from <- to <- pearson <- NULL

    logger::log_success('Running ligand activity analysis.')

    ligand_activites_rds_path <-
        nichenet_results_dir() %>%
        file.path('ligand_activities.rds')

    potential_ligands <-
        lr_network %>%
        filter(
            from %in% expressed_genes_transmitter &
            to %in% expressed_genes_receiver
        ) %>%
        pull(from) %>%
        unique

    genes_of_interest %<>%
        intersect(rownames(ligand_target_matrix)) %>%
        setdiff(potential_ligands)

    background_genes %<>%
        if_null(
            expressed_genes_receiver %>%
            intersect(rownames(ligand_target_matrix))
            # shouldn't we also remove the genes of interest?
        )

    nichenetr::predict_ligand_activities(
        geneset = genes_of_interest,
        background_expressed_genes = background_genes,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
    ) %>%
    arrange(-pearson) %>%
    list(
        ligand_activities = .,
        ligand_target_links = nichenet_ligand_target_links(.)
    ) %T>%
    {saveRDS(.$ligand_activities, ligand_activites_rds_path)} %T>%
    {logger::log_success(
        'Finished running ligand activity analysis, saved to `%s`.',
        ligand_activites_rds_path
    )}

}


#' Compiles a table with weighted ligand-target links
#'
#' A wrapper around \code{nichenetr::get_weighted_ligand_target_links} to
#' compile a data frame with weighted links from the top ligands to their
#' top targets.
#'
#' @param ligand_activities Ligand activity table as produced by
#'     \code{nichenetr::predict_ligand_activities}.
#' @param ligand_target_matrix Ligand-target matrix as produced by
#'     \code{nichenetr::construct_ligand_target_matrix} or the wrapper
#'     around it in the current package:
#'     \code{\link{nichenet_ligand_target_matrix}}.
#' @param genes_of_interest Character vector with the gene symbols of the
#'     genes of interest. These are the genes in the receiver cell population
#'     that are potentially affected by ligands expressed by interacting
#'     cells (e.g. genes differentially expressed upon cell-cell interaction).
#' @param n_top_ligands How many of the top ligands to include in the
#'     ligand-target table.
#' @param n_top_targets For each ligand, how many of the top targets to
#'     include in the ligand-target table.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' expression <- nichenet_expression_data()
#' optimization_results <- nichenet_optimization(networks, expression)
#' nichenet_model <- nichenet_build_model(optimization_results, networks)
#' lt_matrix <- nichenet_ligand_target_matrix(
#'     nichenet_model$weighted_networks,
#'     networks$lr_network,
#'     nichenet_model$optimized_parameters
#' )
#' ligand_activities <- nichenet_ligand_activities(
#'     ligand_target_matrix = lt_matrix,
#'     lr_network = networks$lr_network,
#'     # the rest of the parameters should come
#'     # from your transcriptomics data:
#'     expressed_genes_transmitter = expressed_genes_transmitter,
#'     expressed_genes_receiver = expressed_genes_receiver,
#'     genes_of_interest = genes_of_interest
#' )
#' lt_links <- nichenet_ligand_target_links(
#'     ligand_activities = ligand_activities,
#'     ligand_target_matrix = lt_matrix,
#'     genes_of_interest = genes_of_interest,
#'     n_top_ligands = 20,
#'     n_top_targets = 100
#' )
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr top_n arrange
#' @importFrom purrr map
nichenet_ligand_target_links <- function(
    ligand_activities,
    ligand_target_matrix,
    genes_of_interest,
    n_top_ligands = 42,
    n_top_targets = 250
){

    # NSE vs. R CMD check workaround
    pearson <- NULL

    ligand_activities %>%
    arrange(-pearson)
    top_n(pearson, n_top_ligands) %>%
    map(
        nichenetr::get_weighted_ligand_target_links,
        geneset = genes_of_interest,
        ligand_target_matrix = ligand_target_matrix,
        n = n_top_targets
    )

}


#' Path to the current NicheNet results directory
#'
#' Path to the directory to save intermediate and final outputs from NicheNet
#' methods.
#'
#' @examples
#' \donttest{
#' nichenet_results_dir()
#' # [1] "nichenet_results"
#' }
#'
#' @export
#' @importFrom magrittr %>% %T>%
nichenet_results_dir <- function(){

    'omnipath.nichenet_results_dir' %>%
    options %>%
    `[[`(1) %T>%
    dir.create(showWarnings = FALSE, recursive = TRUE)

}


#' Builds NicheNet network prior knowledge
#'
#' Builds network knowledge required by NicheNet. For this it calls
#' a multitude of methods to download and combine data from various
#' databases according to the settings. The content of the prior knowledge
#' data is highly customizable, see the documentation of the related
#' functions.
#'
#' @param signaling_network A list of parameters for building the signaling
#'     network, passed to \code{\link{nichenet_signaling_network}}
#' @param lr_network A list of parameters for building the ligand-receptor
#'     network, passed to \code{\link{nichenet_lr_network}}
#' @param gr_network A list of parameters for building the gene regulatory
#'     network, passed to \code{\link{nichenet_gr_network}}
#' @param only_omnipath Logical: a shortcut to use only OmniPath as network
#'     resource.
#'
#' @examples
#' \donttest{
#' networks <- nichenet_networks()
#' networks$gr_network %>% dplyr::sample_n(10)
#' # # A tibble: 10 x 4
#' #    from    to       source               database
#' #    <chr>   <chr>    <chr>                <chr>
#' #  1 MAX     ALG3     harmonizome_ENCODE   harmonizome
#' #  2 MAX     IMPDH1   harmonizome_ENCODE   harmonizome
#' #  3 SMAD5   LCP1     Remap_5              Remap
#' #  4 HNF4A   TNFRSF19 harmonizome_CHEA     harmonizome
#' #  5 SMC3    FAP      harmonizome_ENCODE   harmonizome
#' #  6 E2F6    HIST1H1B harmonizome_ENCODE   harmonizome
#' #  7 TFAP2C  MAT2B    harmonizome_ENCODE   harmonizome
#' #  8 USF1    TBX4     harmonizome_TRANSFAC harmonizome
#' #  9 MIR133B FETUB    harmonizome_TRANSFAC harmonizome
#' # 10 SP4     HNRNPH2  harmonizome_ENCODE   harmonizome
#'
#' # use only OmniPath:
#' omnipath_networks <- nichenet_networks(only_omnipath = TRUE)
#' }
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom purrr map2
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_signaling_network}}}
#'     \item{\code{\link{nichenet_lr_network}}}
#'     \item{\code{\link{nichenet_gr_network}}}
#' }
nichenet_networks <- function(
    signaling_network = list(),
    lr_network = list(),
    gr_network = list(),
    only_omnipath = FALSE
){

    environment() %>%
    as.list() %T>%
    {logger::log_success('Building NicheNet network knowledge')} %>%
    discard(names(.) == 'only_omnipath') %>%
    map2(
        names(.),
        function(args, network_type){
            if(!('only_omnipath' %in% names(args))){
                args$only_omnipath <- only_omnipath
            }
            network_type %>%
            sprintf('nichenet_%s', .) %>%
            get() %>%
            do.call(args)
        }
    ) %T>%
    {logger::log_success('Finished building NicheNet network knowledge')}

}


#' Builds a NicheNet signaling network
#'
#' Builds signaling network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_omnipath}}.
#' @param pathwaycommons List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_pathwaycommons}}.
#' @param harmonizome List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_harmonizome}}.
#' @param vinayagam List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_vinayagam}}.
#' @param cpdb List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_cpdb}}.
#' @param evex List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_evex}}.
#' @param inbiomap List with paramaters to be passed to
#'     \code{\link{nichenet_signaling_network_inbiomap}}.
#' @param only_omnipath Logical: a shortcut to use only OmniPath as network
#'     resource.
#'
#' @examples
#' \donttest{
#' # load everything with the default parameters:
#' sig_network <- nichenet_signaling_network()
#'
#' # override parameters for some resources:
#' sig_network <- nichenet_signaling_network(
#'     omnipath = list(resources = c('SIGNOR', 'SignaLink3', 'SPIKE')),
#'     patwaycommons = NULL,
#'     harmonizome = list(datasets = c('phosphositeplus', 'depod')),
#'     cpdb = list(complex_max_size = 1, min_score = .98),
#'     evex = list(min_confidence = 1.5)
#' )
#'
#' # use only OmniPath:
#' sig_network_omnipath <- nichenet_signaling_network(only_omnipath = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_signaling_network_omnipath}}}
#'     \item{\code{\link{nichenet_signaling_network_pathwaycommons}}}
#'     \item{\code{\link{nichenet_signaling_network_harmonizome}}}
#'     \item{\code{\link{nichenet_signaling_network_vinayagam}}}
#'     \item{\code{\link{nichenet_signaling_network_cpdb}}}
#'     \item{\code{\link{nichenet_signaling_network_evex}}}
#'     \item{\code{\link{nichenet_signaling_network_inbiomap}}}
#' }
nichenet_signaling_network <- function(
    omnipath = list(),
    pathwaycommons = list(),
    harmonizome = list(),
    vinayagam = list(),
    cpdb = list(),
    evex = list(),
    inbiomap = list(),
    only_omnipath = FALSE
){

    environment() %>%
    as.list() %>%
    `[[<-`('network_type', 'signaling') %>%
    do.call(nichenet_network, .)

}


#' Builds a NicheNet ligand-receptor network
#'
#' Builds ligand-receptor network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_omnipath}}.
#' @param guide2pharma List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_guide2pharma}}.
#' @param ramilowski List with paramaters to be passed to
#'     \code{\link{nichenet_lr_network_ramilowski}}.
#' @param only_omnipath Logical: a shortcut to use only OmniPath as network
#'     resource.
#'
#' @examples
#' \donttest{
#' # load everything with the default parameters:
#' lr_network <- nichenet_lr_network()
#'
#' # don't use Ramilowski:
#' lr_network <- nichenet_lr_network(ramilowski = NULL)
#'
#' # use only OmniPath:
#' lr_network_omnipath <- nichenet_lr_network(only_omnipath = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_lr_network_omnipath}}}
#'     \item{\code{\link{nichenet_lr_network_guide2pharma}}}
#'     \item{\code{\link{nichenet_lr_network_ramilowski}}}
#' }
nichenet_lr_network <- function(
    omnipath = list(),
    guide2pharma = list(),
    ramilowski = list(),
    only_omnipath = FALSE
){

    environment() %>%
    as.list() %>%
    `[[<-`('network_type', 'lr') %>%
    do.call(nichenet_network, .)

}


#' Builds a NicheNet gene regulatory network
#'
#' Builds gene regulatory network prior knowledge for NicheNet using multiple
#' resources.
#'
#' @param omnipath List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_omnipath}}.
#' @param harmonizome List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_harmonizome}}.
#' @param regnetwork List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_regnetwork}}.
#' @param htridb List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_htridb}}.
#' @param remap List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_remap}}.
#' @param evex List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_evex}}.
#' @param pathwaycommons List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_pathwaycommons}}.
#' @param trrust List with paramaters to be passed to
#'     \code{\link{nichenet_gr_network_trrust}}.
#' @param only_omnipath Logical: a shortcut to use only OmniPath as network
#'     resource.
#'
#' @examples
#' \donttest{
#' # load everything with the default parameters:
#' gr_network <- nichenet_gr_network()
#'
#' # less targets from ReMap, not using RegNetwork:
#' gr_network <- nichenet_gr_network(
#'     remap = list(top_targets = 200),
#'     regnetwork = NULL,
#' )
#'
#' # use only OmniPath:
#' gr_network_omnipath <- nichenet_gr_network(only_omnipath = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_gr_network_evex}}}
# #'     \item{\code{\link{nichenet_gr_network_harmonizome}}}
#'     \item{\code{\link{nichenet_gr_network_htridb}}}
#'     \item{\code{\link{nichenet_gr_network_omnipath}}}
#'     \item{\code{\link{nichenet_gr_network_pathwaycommons}}}
#'     \item{\code{\link{nichenet_gr_network_regnetwork}}}
#'     \item{\code{\link{nichenet_gr_network_remap}}}
#'     \item{\code{\link{nichenet_gr_network_trrust}}}
#' }
nichenet_gr_network <- function(
    omnipath = list(),
    harmonizome = list(),
    regnetwork = list(),
    htridb = list(),
    remap = list(),
    evex = list(),
    pathwaycommons = list(),
    trrust = list(),
    only_omnipath = FALSE
){

    environment() %>%
    as.list() %>%
    `[[<-`('network_type', 'gr') %>%
    do.call(nichenet_network, .)

}


#' Common method to build NicheNet network prior knowledge
#'
#' @param network_type Character: type of the interactions, either
#'     "signaling", "lr" (ligand-receptor) or "gr" (gene regulatory).
#' @param only_omnipath Logical: a shortcut to use only OmniPath as network
#'     resource.
#' @param ... Argument names are the name of the resources to download (all
#'     lowercase), while their values are lists of arguments to the resource
#'     specific nichenet import methods (an empty list if no arguments should
#'     be overridden). If the value is NULL the resource will be omitted.
#'
#' @return A data frame with interactions suitable for use with NicheNet.
#'
#' @examples
#' \donttest{
#' # load the ligand-receptor network with the default parameters:
#' lr_network <- nichenet_network(network_type = 'lr')
#' }
#'
#' @importFrom purrr map2 discard keep
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows filter
#' @importFrom tibble as_tibble
nichenet_network <- function(network_type, only_omnipath = FALSE, ...){

    # NSE vs. R CMD check workaround
    from <- to <- NULL

    network_types <- list(
        signaling = 'signaling',
        lr = 'ligand-receptor',
        gr = 'gene regulatory'
    )

    list(...) %>%
    {`if`(
        length(.) == 0,
        sprintf('nichenet_%s_network', network_type) %>%
        get() %>%
        formals() %>%
        discard(names(.) == 'only_omnipath'),
        .
    )} %>%
    discard(is.null) %>%
    {`if`(
        only_omnipath,
        keep(., names(.) == 'omnipath'),
        .
    )} %T>%
    {logger::log_success(
        'Starting to build NicheNet %s network',
        network_types[[network_type]]
    )} %>%
    map2(
        names(.),
        function(args, resource){
            resource %>%
            sprintf('nichenet_%s_network_%s', network_type, .) %>%
            get() %>%
            do.call(args)
        }
    ) %>%
    # add an empty tibble just to avoid error in case of using no resources
    c(
        list(
            character() %>%
            list() %>%
            rep(4) %>%
            setNames(c('from', 'to', 'source', 'database')) %>%
            do.call(what = tibble)
        )
    ) %>%
    bind_rows %>%
    filter(from != to) %>%
    filter(!is.na(from) & !is.na(to)) %>%
    as_tibble %T>%
    {logger::log_success(
        'Finished building NicheNet %s network: %d records total',
        network_types[[network_type]],
        nrow(.)
    )}

}


#' Builds signaling network for NicheNet using OmniPath
#'
#' Retrieves network prior knowledge from OmniPath and provides it in
#' a format suitable for NicheNet.
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_post_translational_interactions}}
#'
#' @examples
#' \donttest{
#' # use interactions with at least 2 evidences (reference or database)
#' op_signaling_network <- nichenet_signaling_network_omnipath(
#'     min_curation_effort = 2
#' )
#' }
#'
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_signaling_network}}}
#' }
nichenet_signaling_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    args <- list(...)
    args$exclude %<>% union('ligrecextra')
    args$entity_types <- 'protein'

    do.call(import_post_translational_interactions, args) %>%
    omnipath_interactions_postprocess(type = 'signaling')

}


#' Builds ligand-receptor network for NicheNet using OmniPath
#'
#' Retrieves network prior knowledge from OmniPath and provides it in
#' a format suitable for NicheNet.
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_intercell_network}}
#'
#' @examples
#' \donttest{
#' # use only ligand-receptor interactions (not for example ECM-adhesion):
#' op_lr_network <- nichenet_lr_network_omnipath(ligand_receptor = TRUE)
#'
#' # use only CellPhoneDB and Guide to Pharmacology:
#' op_lr_network <- nichenet_lr_network_omnipath(
#'     resources = c('CellPhoneDB', 'Guide2Pharma')
#' )
#'
#' # only interactions where the receiver is a transporter:
#' op_lr_network <- nichenet_lr_network_omnipath(
#'     receiver_param = list(parent = 'transporter')
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_lr_network}}}
#'     \item{\code{\link{import_intercell_network}}}
#' }
nichenet_lr_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    import_intercell_network(...) %>%
    omnipath_interactions_postprocess(type = 'lr')

}


#' Builds gene regulatory network for NicheNet using OmniPath
#'
#' Retrieves network prior knowledge from OmniPath and provides it in
#' a format suitable for NicheNet.
#' This method never downloads the `ligrecextra` dataset because the
#' ligand-receptor interactions are supposed to come from \code{
#' \link{nichenet_lr_network_omnipath}}.
#'
#' @param min_curation_effort Lower threshold for curation effort
#' @param ... Passed to \code{\link{import_transcriptional_interactions}}
#'
#' @examples
#' \donttest{
#' # use interactions up to confidence level "C" from DoRothEA:
#' op_gr_network <- nichenet_gr_network_omnipath(
#'     dorothea_levels = c('A', 'B', 'C')
#' )
#' }
#'
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @seealso \code{\link{nichenet_gr_network_evex},
#'     \link{nichenet_gr_network_harmonizome},
#'     \link{nichenet_gr_network_htridb},
#'     \link{nichenet_gr_network_omnipath},
#'     \link{nichenet_gr_network_pathwaycommons},
#'     \link{nichenet_gr_network_regnetwork},
#'     \link{nichenet_gr_network_remap},
#'     \link{nichenet_gr_network_trrust}}
nichenet_gr_network_omnipath <- function(
    min_curation_effort = 0,
    ...
){

    args <- list(...)
    args$exclude %<>% union('ligrecextra')
    args$entity_types <- 'protein'

    do.call(import_transcriptional_interactions, args) %>%
    omnipath_interactions_postprocess(type = 'gr')

}


#' Processes OmniPath interactions table into NicheNet format
#'
#' @importFrom dplyr select mutate distinct
#' @importFrom tidyr separate_rows
#' @importFrom magrittr %>%
#'
#' @noRd
omnipath_interactions_postprocess <- function(interactions, type){

    # NSE vs. R CMD check workaround
    from <- to <- NULL

    interactions %>%
    select(from = source_genesymbol, to = target_genesymbol, is_directed) %>%
    # expanding complexes
    separate_rows(from, sep = '_') %>%
    separate_rows(to, sep = '_') %>%
    mutate(
        source = sprintf(
            'omnipath_%sdirected_%s',
            ifelse(is_directed, '', 'un'),
            type
        ),
        database = sprintf('omnipath_%s', type)
    ) %>%
    distinct() %>%
    select(-is_directed)

}


#' NicheNet signaling network from PathwayCommons
#'
#' Builds signaling network prior knowledge for NicheNet using PathwayCommons.
#'
#' @param interaction_types Character vector with PathwayCommons interaction
#'     types. Please refer to the default value and the PathwayCommons
#'     webpage.
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' # use only the "controls-transport-of" interactions:
#' pc_signaling_network <- nichenet_signaling_network_pathwaycommons(
#'     interaction_types = 'controls-transport-of'
#' )
#' }
#'
#' @export
nichenet_signaling_network_pathwaycommons <- function(
    interaction_types = c(
        'catalysis-precedes',
        'controls-phosphorylation-of',
        'controls-state-change-of',
        'controls-transport-of',
        'in-complex-with',
        'interacts-with'
    ),
    ...
){

    nichenet_pathwaycommons_common(
        interaction_types = interaction_types,
        label = 'signaling'
    )

}


#' NicheNet signaling network from Harmonizome
#'
#' Builds signaling network prior knowledge for NicheNet using Harmonizome
#'
#' @param datasets The datasets to use. For possible values please refer to
#'     default value and the Harmonizome webpage.
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' # use only KEA and PhosphoSite:
#' hz_signaling_network <- nichenet_signaling_network_harmonizome(
#'     datasets = c('kea', 'phosphositeplus')
#' )
#' }
#'
#' @export
nichenet_signaling_network_harmonizome <- function(
    datasets = c(
        'phosphositeplus',
        'kea',
        'depod'
    ),
    ...
){

    dataset_names <- list(
        phosphositeplus = 'PhosphoSite',
        kea = 'KEA',
        depod = 'DEPOD'
    )

    harmonizome_nichenet(datasets, dataset_names)

}


#' Combines multiple Harmonizome datasets and converts them to NicheNet format
#'
#' @importFrom dplyr mutate bind_rows
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @seealso \code{\link{harmonizome_download},
#'     \link{harmonizome_nichenet_process}}
#'
#' @noRd
harmonizome_nichenet <- function(datasets, dataset_names){

    datasets %>%
    map(harmonizome_nichenet_process) %>%
    bind_rows() %>%
    mutate(
        source = sprintf('harmonizome_%s', dataset_names[source]),
        database = 'harmonizome'
    )

}


#' Processes a table downloaded from Harmonizome to NicheNet format
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang sym
#' @importFrom dplyr select mutate
#' @importFrom stringr str_split_fixed
#' @seealso \code{\link{harmonizome_download}, \link{harmonizome_nichenet}}
#'
#' @noRd
harmonizome_nichenet_process <- function(dataset){

    # NSE vs. R CMD check workaround
    to <- NULL

    target_desc_col <- c('geotf', 'geokinase', 'geogene')
    to_col <- `if`(
        dataset %in% target_desc_col,
        sym('target_desc'),
        sym('target')
    )
    target_proc <- list(
        geotf = toupper,
        msigdbonc = function(x){
            x %>% str_split_fixed('[._]', 2) %>% `[`(,1)
        }
    )

    dataset %>%
    harmonizome_download() %>%
    select(from = source, to = !!to_col) %>%
    {`if`(
        dataset %in% names(target_proc),
        mutate(., to = target_proc[[dataset]](to)),
        .
    )} %>%
    mutate(source = dataset)

}


#' NicheNet signaling network from Vinayagam
#'
#' Builds signaling network prior knowledge for NicheNet using Vinayagam 2011
#' Supplementary Table S6. Find out more at
#' https://doi.org/10.1126/scisignal.2001699
#'
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' vi_signaling_network <- nichenet_signaling_network_vinayagam()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select rename distinct
#' @export
nichenet_signaling_network_vinayagam <- function(...){

    # NSE vs. R CMD check workaround
    `Input-node Gene Symbol` <- `Output-node Gene Symbol` <- NULL

    vinayagam_download() %>%
    nichenet_common_postprocess(
        from_col = `Input-node Gene Symbol`,
        to_col = `Output-node Gene Symbol`,
        source = 'vinayagam_ppi',
        database = 'vinayagam'
    )

}


#' Builds signaling network prior knowledge for NicheNet using ConsensusPathDB
#' (CPDB)
#'
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' # use some parameters stricter than default:
#' cpdb_signaling_network <- nichenet_signaling_network_cpdb(
#'     complex_max_size = 2,
#'     min_score = .99
#' )
#' }
#'
#' @importFrom dplyr select mutate distinct
#' @importFrom magrittr %>%
#' @export
#' @seealso \code{\link{nichenet_signaling_network},
#'     \link{consensuspathdb_download}}
nichenet_signaling_network_cpdb <- function(...){

    # NSE vs. R CMD check workaround
    in_complex <- genesymbol_a <- genesymbol_b <- from <- to <- NULL

    consensuspathdb_download(...) %>%
    mutate(
        source = sprintf(
            'cpdb_%s',
            ifelse(in_complex, 'complex', 'interaction')
        ),
        database = 'cpdb'
    ) %>%
    rename(from = genesymbol_a, to = genesymbol_b) %>%
    select(from, to, source, database) %>%
    distinct()

}


#' NicheNet signaling network from EVEX
#'
#' Builds signaling network prior knowledge for NicheNet from the EVEX
#' database.
#'
#' @param top_confidence Double, between 0 and 1. Threshold based on the
#'     quantile of the confidence score.
#' @param indirect Logical: whether to include indirect interactions.
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' ev_signaling_network <- nichenet_signaling_network_evex(
#'     top_confidence = .9
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate filter
#' @export
#'
#' @seealso \itemize{
#'     \item{\code{\link{evex_download}}}
#'     \item{\code{\link{nichenet_signaling_network}}}
#' }
nichenet_signaling_network_evex <- function(
    top_confidence = .75,
    indirect = FALSE,
    ...
){

    # NSE vs. R CMD check workaround
    coarse_type <- refined_type <- NULL

    categories <- list(
        Binding = 'binding',
        `Regulation of binding` = 'regulation_binding',
        Regulation_of_phosphorylation = 'phosphorylation'
    )

    evex_download(top_confidence = top_confidence, ...) %>%
    select(
        from = source_genesymbol,
        to = target_genesymbol,
        coarse_type,
        refined_type
    ) %>%
    {`if`(
        indirect,
        .,
        filter(., coarse_type != 'Indirect_regulation')
    )} %>%
    filter(
        !(refined_type %in% c(
            # these belong to transcriptional regulation
            'Regulation of expression',
            'Regulation of transcription',
            'Catalysis of DNA methylation'
        ))
    ) %>%
    mutate(
        source = sprintf(
            'evex_%s',
            ifelse(
                refined_type %in% names(categories),
                categories[refined_type],
                ifelse(
                    startsWith(refined_type, 'Catalysis'),
                    'catalysis',
                    ifelse(
                        coarse_type == 'Regulation',
                        'regulation_other',
                        'binding'
                    )
                )
            )
        ),
        database = 'evex_signaling'
    ) %>%
    select(-coarse_type, -refined_type)


}


#' NicheNet signaling network from InWeb InBioMap
#'
#' Builds signaling network prior knowledge for NicheNet from the InWeb
#' InBioMap database.
#'
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' ib_signaling_network <- nichenet_signaling_network_inbiomap()
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
#'
#' @seealso \code{\link{nichenet_signaling_network}, \link{inbiomap_download}}
nichenet_signaling_network_inbiomap <- function(...){

    # NSE vs. R CMD check workaround
    genesymbol_a <- genesymbol_b <- NULL

    inbiomap_download(...) %>%
    nichenet_common_postprocess(
        from_col = genesymbol_a,
        to_col = genesymbol_b,
        source = 'inweb_interaction',
        database = 'inweb_inbiomap'
    )

}


#' Ligand-receptor network from Guide to Pharmacology
#'
#' Downloads ligand-receptor interactions from the Guide to Pharmacology
#' database and converts it to a format suitable for NicheNet.
#'
#' @return Data frame with ligand-receptor interactions in NicheNet format.
#'
#' @examples
#' \donttest{
#' g2p_lr_network <- nichenet_lr_network_guide2pharma()
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \code{\link{nichenet_lr_network}, \link{guide2pharma_download}}
nichenet_lr_network_guide2pharma <- function(){

    # NSE vs. R CMD check workaround
    target_species <- ligand_species <- ligand_gene_symbol <-
        target_gene_symbol <- NULL

    guide2pharma_download() %>%
    filter(
        target_species == 'Human' &
        ligand_species == 'Human'
    ) %>%
    nichenet_common_postprocess(
        source = 'pharmacology',
        database = 'guide2pharmacology',
        from_col = ligand_gene_symbol,
        to_col = target_gene_symbol
    )

}


#' Ligand-receptor network from Ramilowski 2015
#'
#' Downloads ligand-receptor interactions from Supplementary Table 2 of the
#' paper 'A draft network of ligand–receptor-mediated multicellular signalling
#' in human' (Ramilowski et al. 2015,
#' https://www.nature.com/articles/ncomms8866). It converts the downloaded
#' table to a format suitable for NicheNet.
#'
#' @param evidences Character: evidence types, "literature supported",
#'     "putative" or both.
#'
#' @return Data frame with ligand-receptor interactions in NicheNet format.
#'
#' @examples
#' \donttest{
#' # use only the literature supported data:
#' rami_lr_network <- nichenet_lr_network_ramilowski(
#'     evidences = 'literature supported'
#' )
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_lr_network}}}
#'     \item{\code{\link{ramilowski_download}}}
#' }
nichenet_lr_network_ramilowski <- function(
    evidences = c('literature supported', 'putative')
){

    # NSE vs. R CMD check workaround
    Pair.Evidence <- Ligand.ApprovedSymbol <- Receptor.ApprovedSymbol <- NULL

    ramilowski_download() %>%
    filter(Pair.Evidence %in% evidences) %>%
    nichenet_common_postprocess(
        source = 'ramilowski_known',
        database = 'ramilowski',
        from_col = Ligand.ApprovedSymbol,
        to_col = Receptor.ApprovedSymbol
    )

}


#' NicheNet gene regulatory network from Harmonizome
#'
#' Builds gene regulatory network prior knowledge for NicheNet using
#' Harmonizome
#'
#' @param datasets The datasets to use. For possible values please refer to
#'     default value and the Harmonizome webpage.
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' # use only JASPAR and TRANSFAC:
#' hz_gr_network <- nichenet_gr_network_harmonizome(
#'     datasets = c('jasparpwm', 'transfac', 'transfacpwm')
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_gr_network}}}
#'     \item{\code{\link{harmonizome_download}}}
#' }
nichenet_gr_network_harmonizome <- function(
    datasets = c(
        'cheappi',
        'encodetfppi',
        'jasparpwm',
        'transfac',
        'transfacpwm',
        'motifmap',
        'geotf',
        'geokinase',
        'geogene'
    ),
    ...
){

    # NSE vs. R CMD check workaround
    to <- from <- NULL

    dataset_names <- list(
        cheappi = 'CHEA',
        encodetfppi = 'ENCODE',
        jaspar = 'JASPAR',
        transfac = 'TRANSFAC_CUR',
        transfacpwm = 'TRANSFAC',
        motifmap = 'MOTIFMAP',
        geotf = 'GEO_TF',
        geokinase = 'GEO_KINASE',
        geogene = 'GEO_GENE',
        msigdbonc = 'MSIGDB_GENE'
    )

    harmonizome_nichenet(datasets, dataset_names) %>%
    rename(from = to, to = from)

}


#' NicheNet gene regulatory network from RegNetwork
#'
#' Builds a gene regulatory network using data from the RegNetwork database
#' and converts it to a format suitable for NicheNet.
#'
#' @examples
#' \donttest{
#' regn_gr_network <- nichenet_gr_network_regnetwork()
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @seealso \itemize{
#'     \item{\code{\link{regnetwork_download}}}
#'     \item{\code{\link{nichenet_gr_network}}}
#' }
nichenet_gr_network_regnetwork <- function(){

    # NSE vs. R CMD check workaround
    source_type <- target_type <- NULL

    regnetwork_download() %>%
    filter(
        source_type == 'protein' &
        target_type == 'protein'
    ) %>%
    nichenet_common_postprocess(
        source = 'regnetwork_source',
        database = 'regnetwork'
    )

}


#' NicheNet gene regulatory network from TRRUST
#'
#' Builds a gene regulatory network using data from the TRRUST database
#' and converts it to a format suitable for NicheNet.
#'
#' @examples
#' \donttest{
#' trrust_gr_network <- nichenet_gr_network_trrust()
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \itemize{
#'     \item{\code{\link{trrust_download}}}
#'     \item{\code{\link{nichenet_gr_network}}}
#' }
nichenet_gr_network_trrust <- function(){

    trrust_download() %>%
    nichenet_common_postprocess(
        source = 'trrust',
        database = 'trrust'
    )

}


#' NicheNet gene regulatory network from HTRIdb
#'
#' Builds a gene regulatory network using data from the HTRIdb database
#' and converts it to a format suitable for NicheNet.
#'
#' @examples
#' \donttest{
#' htri_gr_network <- nichenet_gr_network_htridb()
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{htridb_download}, \link{nichenet_gr_network}}
nichenet_gr_network_htridb <- function(){

    # NSE vs. R CMD check workaround
    SYMBOL_TF <- SYMBOL_TG <- NULL

    htridb_download() %>%
    nichenet_common_postprocess(
        source = 'HTRIDB',
        database = 'HTRIDB',
        from_col = SYMBOL_TF,
        to_col = SYMBOL_TG
    )

}


#' NicheNet gene regulatory network from ReMap
#'
#' Builds a gene regulatory network using data from the ReMap database
#' and converts it to a format suitable for NicheNet.
#'
#' @param score Numeric: a minimum score between 0 and 1000, records with
#'     lower scores will be excluded. If NULL no filtering performed.
#' @param top_targets Numeric: the number of top scoring targets for each
#'     TF. Essentially the maximum number of targets per TF. If NULL the
#'     number of targets is not restricted.
#' @param only_known_tfs Logical: whether to exclude TFs which are not in
#'     TF census.
#'
#' @examples
#' \donttest{
#' # use only max. top 100 targets for each TF:
#' remap_gr_network <- nichenet_gr_network_remap(top_targets = 100)
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \itemize{
#'     \item{\code{\link{remap_filtered}}}
#'     \item{\code{\link{nichenet_gr_network}}}
#' }
nichenet_gr_network_remap <- function(
    score = 100,
    top_targets = 500,
    only_known_tfs = TRUE
){

    remap_filtered(
        score = score,
        top_targets = top_targets,
        only_known_tfs = only_known_tfs
    ) %>%
    nichenet_common_postprocess(
        source = 'Remap_5',
        database = 'Remap'
    )

}


#' NicheNet gene regulatory network from EVEX
#'
#' Builds a gene regulatory network using data from the EVEX database
#' and converts it to a format suitable for NicheNet.
#'
#' @return Data frame of interactions in NicheNet format.
#'
#' @param top_confidence Double, between 0 and 1. Threshold based on the
#' quantile of the confidence score.
#' @param indirect Logical: whether to include indirect interactions.
#' @param regulation_of_expression Logical: whether to include also the
#'     "regulation of expression" type interactions.
#'
#' @examples
#' \donttest{
#' # use only the 10% with the highest confidence:
#' evex_gr_network <- nichenet_gr_network_evex(top_confidence = .9)
#' }
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom dplyr filter select mutate
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_gr_network}}}
#'     \item{\code{\link{evex_download}}}
#' }
nichenet_gr_network_evex <- function(
    top_confidence = .75,
    indirect = FALSE,
    regulation_of_expression = FALSE
){

    # NSE vs. R CMD check workaround
    confidence <- coarse_type <- refined_type <- NULL

    gr_types <- `if`(
        regulation_of_expression,
        c('Regulation of expression', 'Regulation of transcription'),
        'Regulation of transcription'
    )

    evex_download() %>%
    filter(confidence > quantile(confidence, top_confidence)) %>%
    {`if`(
        indirect,
        .,
        filter(., coarse_type != 'Indirect_regulation')
    )} %>%
    filter(
        refined_type %in% gr_types
    ) %>%
    nichenet_common_postprocess(
        source = 'evex_regulation_expression',
        database = 'evex_expression'
    )

}


#' NicheNet gene regulatory network from PathwayCommons
#'
#' Builds gene regulation prior knowledge for NicheNet using PathwayCommons.
#'
#' @param interaction_types Character vector with PathwayCommons interaction
#'     types. Please refer to the default value and the PathwayCommons
#'     webpage.
#' @param ... Ignored.
#'
#' @examples
#' \donttest{
#' pc_gr_network <- nichenet_gr_network_pathwaycommons()
#' }
#'
#' @export
#' @seealso \itemize{
#'     \item{\code{\link{nichenet_gr_network}}}
#'     \item{\code{\link{pathwaycommons_download}}}
#' }
nichenet_gr_network_pathwaycommons <- function(
    interaction_types = 'controls-expression-of',
    ...
){

    nichenet_pathwaycommons_common(
        interaction_types = interaction_types,
        label = 'expression'
    )

}


#' Retrieves interactions from PathwayCommons and converts them to NicheNet
#' format
#'
#' @param interaction_types Character vector with PathwayCommons interaction
#'     types. Please refer to the default value and the PathwayCommons
#'     webpage.
#' @param label Character: suffix for the NicheNet `database` field:
#'     "signaling" for the signaling network and "expression" for gene
#'     regulatory network.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter relocate select
#'
#' @noRd
nichenet_pathwaycommons_common <- function(interaction_types, label){

    # NSE vs. R CMD check workaround
    type <- from <- to <- NULL

    pathwaycommons_download() %>%
    filter(
        type %in% interaction_types
    ) %>%
    mutate(
        source = sprintf(
            'pathwaycommons_%s',
            gsub('-', '_', type, fixed = TRUE)
        ),
        database = sprintf('pathwaycommons_%s', label)
    ) %>%
    relocate(from, to) %>%
    select(-type)

}


#' Common postprocessing from building a NicheNet format network table
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !! enquo
#' @importFrom dplyr select distinct mutate
#'
#' @noRd
nichenet_common_postprocess <- function(
    data,
    source,
    database,
    from_col = source_genesymbol,
    to_col = target_genesymbol
){

    from_col <- enquo(from_col)
    to_col <- enquo(to_col)

    data %>%
    select(
        from = !!from_col,
        to = !!to_col
    ) %>%
    distinct() %>%
    mutate(
        source = source,
        database = database
    )

}


#' Expression data from ligand-receptor perturbation experiments used by
#' NicheNet
#'
#' NicheNet uses expression data from a collection of published ligand or
#' receptor KO or perturbation experiments to build its model. This function
#' retrieves the original expression data, deposited in Zenodo
#' (https://zenodo.org/record/3260758).
#'
#' @return Nested list, each element contains a data frame of processed
#'     expression data and key variables about the experiment.
#'
#' @examples
#' \donttest{
#' exp_data <- nichenet_expression_data()
#' names(exp_data) %>% head()
#' # [1] "bmp4_tgfb"     "tgfb_bmp4"     "nodal_Nodal"   "spectrum_Il4"
#' # [5] "spectrum_Tnf"  "spectrum_Ifng"
#' exp_data %>% head() %>% purrr::map_chr('from')
#' #     bmp4_tgfb     tgfb_bmp4   nodal_Nodal  spectrum_Il4  spectrum_Tnf
#' #       "BMP4"       "TGFB1"       "NODAL"         "IL4"         "TNF"
#' # spectrum_Ifng
#' #       "IFNG"
#' }
#'
#' @importFrom magrittr %T>%
#' @export
nichenet_expression_data <- function(){

    generic_downloader(
        url_key = 'omnipath.nichenet_expression_url',
        reader = url_rds,
        reader_param = list(),
        resource = 'NicheNet expression data'
    ) %T>%
    load_success()

}