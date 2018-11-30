#' Compare prediction CIs to actual/expected values; label prediction
#' as accurate if actual value is > 0 and CI does not contain 0 or if
#' actual value == 0 and CI does include 0
#'
#'
bin_predict <- function(
  model, dtf, response_var = extract_rv(model),
  threshold_method = c('CI', 'optimized')) {

  threshold_method <- match.arg(threshold_method,
    choices = c('CI', 'optimized'), several.ok = FALSE)

  stopifnot(response_var %in% colnames(dtf))

  knn_obj <- extract_fit_parsnip(model)$fit

  ## 'raw_CI' prediction is so far only available in my fork of the
  ## kknn package (github.com/slagtermaarten/kknn)
  preds <- predict(knn_obj, dtf, type = 'raw_CI') %>%
    dplyr::bind_cols(dplyr::select(dtf, condition_name)) %>%
    simplify_condition_name() %>%
    dplyr::mutate(
      actual = norm_regressor(vn = response_var, dtf[[response_var]])
    ) %>%
    dplyr::mutate(CI_include_zero = CI_l <= 0 & CI_h < 1) %>%
    dplyr::mutate(CI_include_one = CI_l > 0 & CI_h >= 1) %>%
    dplyr::mutate(CI_non_zero = CI_l > 0) %>%
    dplyr::mutate(CI_zero = CI_h <= 0) %>%
    dplyr::mutate(
      CI_includes_actual = CI_l <= actual & CI_h >= actual) %>%
    dplyr::mutate(
      actual_bin = factor(actual > 0,
        levels = c(FALSE, TRUE), labels = c('0', '> 0'))
    ) %>%
    dplyr::mutate(
      CI_include_zero = factor(CI_include_zero,
        levels = c(TRUE, FALSE),
        labels = c('CI includes 0', 'CI does not include 0'))
    ) %>%
    dplyr::mutate(
      CI_includes_actual = factor(CI_includes_actual,
        levels = c(TRUE, FALSE),
        labels = c('CI spans actual', 'CI does not span actual'))
    ) %>%
    { . }

  ## Annotate the 'accurate' column, indicating whether the particular
  ## cell/sample in that particular row is on the correct side of the
  ## binary classification threshold
  if (threshold_method == 'optimized')
    preds <- optimize_bin_threshold(preds)
  else if (threshold_method == 'CI')
    preds <- annotate_bin_CI_desirability(preds)

  if (F)
    with(preds, mean(accurate))

  attr(preds, 'response_var') <- response_var

  return(preds)
}


#' Compute 'binarized' accuracy. What percentage of the cases that are
#' supposed to be 0 (>0) are indeed robustly 0 (>0)?
#'
#' DEPRECATED
bin_CI_accuracy <- function(preds, group_vars = 'condition_name') {
  freq_table <- preds %>%
    dplyr::group_by(condition_name, actual_bin, CI_include_zero) %>%
    dplyr::summarize(N = n()) %>%
    dplyr::mutate(freq = N / sum(N)) %>%
    { . }

  freq_table %>%
    annotate_CI_desirability() %>%
    dplyr::group_by_at(group_vars) %>%
    dplyr::summarize(bin_accuracy =
      sum(rep_zero(N[accurate == TRUE])) / sum(N)) %>%
    dplyr::mutate(label = glue::glue('% binarily accurate: \\
      {scales::percent(bin_accuracy)}')) %>%
    dplyr::mutate(label = as.character(label))
}


#' Compute 'binarized' accuracy. What percentage of the cases that are
#' supposed to be 0 (>0) are indeed 0 (>0)? Assess this using a
#' (adapted) threshold, stored in column 'optim_thresh'.
#'
#'
bin_accuracy <- function(preds, group_vars = 'condition_name') {
  freq_table <- preds %>%
    dplyr::group_by(condition_name, actual_bin, y_bin) %>%
    dplyr::summarize(N = n()) %>%
    dplyr::mutate(freq = N / sum(N)) %>%
    dplyr::mutate(accurate = y_bin == (actual_bin == '> 0')) %>%
    { . }

  freq_table %>%
    dplyr::group_by_at(group_vars) %>%
    dplyr::summarize(bin_accuracy =
      sum(rep_zero(N[accurate == TRUE])) / sum(N)) %>%
    dplyr::mutate(label = glue::glue('% binarily accurate: \\
      {scales::percent(bin_accuracy)}')) %>%
    dplyr::mutate(label = as.character(label))
}


mae_accuracy <- function(preds, group_vars = 'condition_name') {
  mae_table <- preds %>%
    dplyr::group_by_at(group_vars) %>%
    dplyr::summarize(mae = yardstick::mae(across(), actual, y)) %>%
    unnest(mae) %>%
    dplyr::select(condition_name, mae = .estimate) %>%
    dplyr::mutate(label = glue::glue('MAE: {round(mae, 2)}')) %>%
    { . }
  return(mae_table)
}


#' Plot scores for randomly selected 'illustrative' cells
#'
#'
plot_binary_score_mod <- function(
  preds,
  response_var = attr(preds, 'response_var'),
  mod_idx = 1,
  fn_app = NULL) {

  ## The desired scores for each of the detected groups
  exp_values <- preds %>%
    dplyr::distinct(condition_name, actual)

  stopifnot(!is.null(response_var))
  nrv <- glue::glue('{response_var}_norm')
  CI_part <- if (any(!is.na(preds$CI_l))) ' with 95% CI' else ''

  p <- preds %>%
    dplyr::group_by(condition_name) %>%
    dplyr::sample_n(100L) %>%
    dplyr::arrange(y) %>%
    dplyr::mutate(idx = 1:n()) %>%
    ggplot(
      mapping = aes(x = idx, y = y, ymin = CI_l, ymax = CI_h,
        colour = accurate)) +
    geom_hline(
      yintercept = preds$optim_thresh[1],
      inherits.aes = F, color = 'blue', linetype = 'dashed'
    ) +
    geom_hline(
      data = exp_values, mapping = aes(yintercept = actual),
      inherits.aes = F, color = 'indianred3'
    ) +
    geom_pointrange(
      position = 'dodge', alpha = .5,
      size = .5, fatten = 1
    ) +
    facet_wrap(~condition_name) +
    scale_x_discrete(name = '', breaks = c()) +
    scale_y_continuous(
      name = glue::glue('Predicted {axis_labels[nrv]}{CI_part}'),
      limits = c(0, 1)
    ) +
    # scale_colour_discrete(name = 'Actual (normalized) value')
    scale_colour_discrete(name = 'Binarily accurate')

  bin_dtf <-
    bin_accuracy(preds, 'condition_name') %>%
    dplyr::filter(!is.na(bin_accuracy))
  if (nrow(bin_dtf) > 0)
    p <- p + ggpp::geom_text_npc(
      data = bin_dtf,
      size = 2, hjust = 0, vjust = 1,
      mapping = aes(npcx = .05, npcy = .95, label = label)
    )

  mae_dtf <-
    mae_accuracy(preds, 'condition_name') %>%
    dplyr::filter(!is.na(mae))
  if (nrow(mae_dtf) > 0)
    p <- p + ggpp::geom_text_npc(
      data = mae_dtf,
      size = 2, hjust = 0, vjust = 1,
      mapping = aes(npcx = .05, npcy = .87, label = label)
    )

  mod_idx <- as.character(mod_idx)
  fn <- file.path(Sys.getenv('img_dir'),
    glue::glue('intuitive_prediction_{response_var}\\
      {make_flag(mod_idx)}{make_flag(fn_app)}.png'))
  print_plot_eval(print(p), filename = fn)
  # np_mod <- train_np(bm_dtf, rv = response_var)
  # sc_dtf$pred <- predict(model, sc_dtf, type = 'raw')
  # p <- maartenutils::plot_scatter_cor(
  #   x_var = as.character(response_var), y_var = 'pred',
  #   dtf = preds)
  # print_plot_eval(print(p),
  #   filename = file.path(Sys.getenv('img_dir'),
  #     glue::glue('actual_vs_pred_{response_var}.png')))

  return(invisible(fn))
}


assess_binary_prediction_accuracy <- function(
  ref_obj,
  query_obj,
  response_vars = 'tnf_duration',
  model_hyperparam_grid = tar_read(model_hyperparam_grid),
  mod_idxs = 1:nrow(model_hyperparam_grid),
  var_regex = '.*',
  N_spec_levels = 6L,
  make_plots = FALSE) {

  if (maartenutils::null_dat(query_obj) || 
      maartenutils::null_dat(ref_obj)) 
    return(NULL)

  if (is.null(model_hyperparam_grid))
    model_hyperparam_grid <-
      gen_knn_hyperparam_grid(N_spec_levels = N_spec_levels)

  if (!missing(response_vars) && !is.null(var_regex)) {
    response_vars <-
      stringr::str_subset(tar_read(response_vars), var_regex)
  }

  response_vars <-
    response_vars %>%
    purrr::keep(~column_present(query_obj, .x))

  if (length(response_vars) == 0)
    return(NULL)

  binary_model_assessment <-
    tidyr::crossing(
      response_var = response_vars,
      mod_idx = mod_idxs
    ) %>%
    dplyr::mutate(
      preds = purrr::map2(response_var, mod_idx,
        function(response_var, mod_idx) {
          out <- tryCatch({
            model <- train_knn_reg_mod(
              dtf = ref_obj,
              response_var = response_var,
              N_spec_levels = N_spec_levels,
              mod_idx = mod_idx
            )
            out <- bin_predict(
              model = model,
              dtf = query_obj,
              threshold_method = 'optimized'
            )
            return(out)
          }, error = function(e) { print(e); NULL })
        })
      )

  if (maartenutils::null_dat(binary_model_assessment))
    return(NULL)

  if (make_plots)
    plot_binary_score_mod(preds, mod_idx = mod_idx)

  binary_model_assessment <-
    binary_model_assessment %>%
    dplyr::mutate(pred_score = purrr::map(preds, bin_accuracy))

  ## Mean per explanatory variable
  mean_assess <-
    binary_model_assessment %>%
    group_by(mod_idx, response_var) %>%
    summarize(acc = map(pred_score, ~.x$bin_accuracy))

  ## Mean of means
  grand_mean <- mean_assess %>%
    dplyr::group_by(mod_idx) %>%
    dplyr::summarize(
      grand_mean = mean(unlist(acc))
    )

  out <- mean_assess %>%
    dplyr::left_join(grand_mean, by = 'mod_idx')

  return(out)
}


optimize_bin_threshold <- function(
  preds,
  group_vars = 'condition_name') {

  #' Where to binarize the predicted score to best separate true zero
  #' and true non-zero values?
  find_optimal_threshold <- function(dtf) {
    if (maartenutils::null_dat(dtf)) return(NA_real_)
    dtf <- dtf %>%
      dplyr::filter(!is.na(y) & !is.na(actual))
    init <- dtf %>%
      dplyr::filter(!is.na(y)) %>%
      dplyr::group_by(actual_bin) %>%
      dplyr::summarize(m = mean(y)) %>%
      dplyr::summarize(m[1] + diff(m) / 2) %>%
      unlist()
    if (length(init) == 0 || !is.finite(init))
      return(tibble(value = unique(mean(dtf$y))))
    dtf$response <- dtf$actual > 0
    if (all(dtf$response == F) || all(dtf$response == T)) {
      return(tibble(value = .5))
    }
    optim_func <- function(th) {
      dtf$y_bin <- dtf$y > th
      # with(dtf, fisher.test(y_bin, response)$p.value)
      with(preds, mean(dtf$y_bin != as.logical(dtf$response)))
    }
    candidates <- seq(0, .9, by = 1e-3)
    scores <- vapply(candidates, optim_func, numeric(1))
    if (maartenutils::eps(var(scores), 0, .Machine$double.eps))
      return(tibble(value = .5))
    out <- broom::tidy(optim(fn = optim_func,
        # par = init,
        par = candidates[dplyr::last(which(scores == min(scores)))],
        lower = 1e-9, upper = 1-1e-9, method = 'L-BFGS-B'))
    # optim_func(out$value)
    # optim_func(candidates[last(which(scores == min(scores)))])
    return(out)
  }
  optim_thresh <- find_optimal_threshold(preds)$value

  preds <- preds %>%
    dplyr::mutate(optim_thresh = optim_thresh) %>%
    dplyr::mutate(y_bin = y > optim_thresh) %>%
    dplyr::mutate(accurate = y_bin == (actual > 0)) %>%
    { . }
  # mean(preds$accurate)

  # optimal_thresholds <-
  #   preds %>%
  #   dplyr::group_by_at(group_vars) %>%
  #   dplyr::summarize(y = mean(y)) %>%
  #   dplyr::mutate(
  #     # summary = summary(data$y),
  #     # optimal_threshold = purrr::map_dbl(data,
  #       # ~nrow(.x)))
  #     # optimal_threshold = purrr::map(data,
  #     optimal_threshold = purrr::map(data,
  #         ~find_optimal_threshold(.x)))
  #     # optimal_threshold = purrr::map_dbl(data,
  #     #     ~nrow(.x)))
  # head(optimal_thresholds)

  return(preds)
}


annotate_bin_CI_desirability <- function(preds) {
  factor_combination_labellings <-
    preds %>%
    dplyr::ungroup() %>%
    dplyr::select(actual_bin, CI_include_zero) %>%
    tidyr::crossing() %>%
    dplyr::mutate(accurate =
      (actual_bin == '0' & CI_include_zero == 'CI includes 0') |
      (actual_bin == '> 0' & CI_include_zero == 'CI does not include 0')
    )

  preds <-
    dplyr::left_join(
      preds,
      factor_combination_labellings,
      by = c('actual_bin', 'CI_include_zero')
    )

  preds$y_bin <- preds$CI_include_zero == 'CI does not include 0'

  attr(preds, 'factor_combination_labellings') <-
    factor_combination_labellings

  return(preds)
}
