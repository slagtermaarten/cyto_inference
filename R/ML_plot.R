cost_labels <- c(
  'prediction_cost' = 'Average',
  'duration_cost' = 'Duration',
  'tnf_conc_cost' = 'TNFa concentration',
  'ifn_conc_cost' = 'IFNy concentration',
  'sn_dilution_cost' = 'T-cell SN dilution'
)


#'
#'
#' @param mod_id_vars Column names pertaining to model type
read_in_benchmarks <- function(
  y_var = 'prediction_cost',
  emb_type = 'ratio',
  mod_id_vars = c('mod_type', 'CV', 'AE'),
  order_var = 'prediction_cost_med',
  experiment = NULL,
  top_n = 3) {
  library(dplyr)

  mod_benchmarks <- controlled_target_read(
    target_name = glue::glue('benchmark_{emb_type}_id_stim')) %>%
    dplyr::mutate(condition_name = droplevels(condition_name))

  rand_benchmarks <- controlled_target_read(
    target_name = glue::glue('benchmark_{emb_type}_random')) %>%
    # dplyr::filter(experiment == .env[['experiment']])
    { . }

  mod_benchmarks %>% tally(experiment, condition_name)
  rand_benchmarks %>% tally(experiment, condition_name)

  ## Among all possible model ID variables, which ones are detected
  ## for this dataset?
  mod_id_vars <- purrr::map_lgl(mod_id_vars,
    ~!is.null(mod_benchmarks[[.x]]) &&
      !all(is.na(mod_benchmarks[[.x]])) &&
    data.table::uniqueN(mod_benchmarks[[.x]]) > 1) %>%
    { mod_id_vars[.] }

  ## Construct model ID name
  mod_benchmarks$mod_id <-
    apply(dplyr::select(mod_benchmarks, any_of(mod_id_vars)), 1,
      function(x) paste(x, collapse = ' '))

  p_dat <-
    dplyr::bind_rows(
      mod_benchmarks,
      rand_benchmarks,
      dplyr::bind_rows(mod_benchmarks, rand_benchmarks) %>%
        group_by(mod_id, mod_type, experiment) %>%
        dplyr::summarize(across(where(is.numeric), mean)) %>%
        dplyr::mutate(condition_name = 'Across conditions')
    ) %>%
    dplyr::mutate(
      mod_id = factor(mod_id,
        levels = naturalsort::naturalsort(unique(as.character(mod_id)))
      )
    )

  allowed_mod_ids <- p_dat %>%
    dplyr::filter(condition_name == 'Across conditions') %>%
    dplyr::group_by(mod_type, experiment) %>%
    dplyr::arrange(prediction_cost_med, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(mod_id) %>%
    as.character() %>%
    { . }
  p_dat <- p_dat %>%
    dplyr::filter(mod_id %in% c('random', allowed_mod_ids)) %>%
    # dplyr::mutate(mod_id = factor(mod_id, 
    #     levels = unique(c('random', allowed_mod_ids))))
    { . }

  p_dat$condition_name <-
    order_factor(p_dat$condition_name,
      exp_levels = levels(mod_benchmarks$condition_name))

  p_dat %>%
    tally(experiment, condition_name, mod_type)

  mod_id_order <- p_dat %>%
    dplyr::filter(condition_name == 'Across conditions') %>%
    dplyr::arrange(.data[[order_var]]) %>%
    dplyr::pull(mod_id) %>%
    as.character() %>%
    setdiff('random') %>%
    unique() %>%
    naturalsort::naturalsort()
  p_dat$mod_id <- order_factor(p_dat$mod_id,
    exp_levels = mod_id_order)

  med_name <- glue::glue('{y_var}_med')
  q9_name <- glue::glue('{y_var}_q9')
  # p_dat <- p_dat %>%
  #   group_by(experiment, condition_name) %>%
  #   dplyr::mutate(
  #     improvement_over_random = .data[[med_name]][mod_id == 'random'])
  # p_dat <- 
  #   q9_name %in% colnames(p_dat)
  # p_dat %>%
  #   group_by(experiment, condition_name) %>%
  #   # head(n = 2) %>%
  #   dplyr::mutate(
  #     improvement_over_random = 
  #       unlist(.data[[med_name]][mod_id == 'random'])[1])

  p_dat <- p_dat %>%
    group_by(experiment, condition_name) %>%
    # head(n = 2) %>%
    dplyr::mutate(
      improvement_over_random = .data[[q9_name]] <
      unlist(.data[[med_name]][mod_id == 'random'])[1])

  p_dat %>% tally(experiment, condition_name, improvement_over_random)

  if (!is.null(experiment)) {
    p_dat <- p_dat %>%
      dplyr::mutate(experiment = as.character(experiment)) %>%
      dplyr::filter(experiment == as.character(.env[['experiment']]))

    p_dat %>% tally(experiment, condition_name, mod_type)
  }

  return(p_dat)
}



#'
#'
#' @param m
plot_benchmark_score_pointranges <- function(
  experiment,
  y_vars = c('prediction_cost', 'duration_cost', 'tnf_conc_cost',
    'ifn_conc_cost', 'sn_dilution_cost')) {

  library(SummarizedExperiment)
  library(ggplot2)

  p_dat <- read_in_benchmarks(top_n = 3, 
    emb_type = 'ratio', experiment = experiment)

  o_fns <- purrr::map(y_vars, function(y_var) {
    med_name <- glue::glue('{y_var}_med')

    if (!med_name %in% colnames(p_dat) ||
        all(is.na(p_dat[[med_name]])))
      return(NULL)

    q1_name <- glue::glue('{y_var}_q1')
    q9_name <- glue::glue('{y_var}_q9')

    p <- ggplot(p_dat,
      aes_string(
        x = 'mod_id', y = med_name,
        # colour = 'improvement_over_random',
        ymin = q1_name, ymax = q9_name)) +
      ylab('Prediction cost') +
      scale_x_discrete('Model ID') +
      # scale_colour_discrete(name = 'Model ID', show.guide = F) +
      geom_pointrange(size = .3, alpha = .8) +
      # stat_summary(fun = 'mean', colour = 'red', geom = 'point') +
      facet_wrap(~condition_name) +
      coord_flip() +
      ggtitle(cost_labels[y_var])

    o_fn <- file.path(exp_plot_dir(experiment),
        glue::glue('prediction_benchmark-\\
          {experiment}\\
          {make_flag(y_var)}.pdf'))
    print_plot_eval(print(p),
      filename = o_fn, width = 12, height = 12, units = 'cm')

    return(o_fn)
  })
}
