limma_preproc <- function(so, down_scale = 1) {
  M <- round(sweep(so2M(so), 1, down_scale, '/'))
  # head(so@meta.data)
  pacman::p_load('limma')
  pacman::p_load('edgeR')
  pacman::p_load('Glimma')
  d0 <- DGEList(counts = M)
  d0 <- calcNormFactors(d0)
  gene_max_e <- apply(edgeR::cpm(d0), 1, max)
  gene_sd <- apply(edgeR::cpm(d0), 1, sd)
  gene_max_e_rc <- apply(M, 1, max)
  # cor(gene_max_e, gene_max_e_rc)
  drop <- which(gene_max_e < 25 | sqrt(gene_sd) < 0.05)
  drop <- which(sqrt(gene_sd) < 0.05)
  d0 <- d0[-drop, ]
  return(d0)
}


estimate_limma <- function(d0, mm, span = .3, extra_fn = '') {
  # diagnose_mm_problems()
  print_plot_eval({ v <- limma::voom(d0, mm, plot=T, span = span) },
    filename = file.path(Sys.getenv('img_dir'),
      glue::glue('limma{extra_fn}.png')))
  v <- limma::voom(d0, mm, plot=F, span = span)
  vfit <- lmFit(v, mm)
  efit <- eBayes(vfit)
  return(efit)
}


limma_model_6369 <- function(so, down_scale = 1, 
  include_interaction = F) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data
  stopifnot(sa$condition_name == colnames(d0))
  if (include_interaction) {
    sa <- separate_duration(sa)
    sa$ifn <- sa$ifn_duration > 0
    sa$ifn_duration <- numeric2factor(sa$ifn_duration)
    sa$ifn_conc <- numeric2factor(sa$ifn_conc)
    sa$duration <- numeric2factor(sa$duration)

    mm <- model.matrix(~1 + ifn_conc + duration + ifn:duration, data = sa)

    ## Leave out the unstimulated sample, such that the intercept of the
    ## model will mostly describe the unstimulated sample
    mm <- mm[, !grepl('^ifn_conc0$', colnames(mm))]
    mm <- mm[, !grepl('^ifn_conc1$', colnames(mm))]
  } else {
    mm <- model.matrix(~1 + ifn_conc + duration, data = sa)
  }

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_6489 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data
  stopifnot(sa$condition_name == colnames(d0))

  if (F) {
    sa <- separate_duration(sa)
    sa$sn <- sa$sn_duration > 0
    sa$sn_duration <- numeric2factor(sa$sn_duration)
    sa$sn_dilution <- numeric2factor(sa$sn_dilution)
    sa$duration <- numeric2factor(sa$duration)
    mm <- model.matrix(
      ~1 + sn_dilution + duration + sn:duration,
      data = sa
    )
    # mm <- mm[, !(colnames(mm) %in% c('sn_dilution0'))]
    mm <- mm[, !(colnames(mm) %in% c('sn_dilution0.0016'))]
    mm <- mm[, !(colnames(mm) %in% c('sn_dilution0.008'))]
    mm <- mm[, !(colnames(mm) %in% c('sn_dilution0.04'))]
    # mm <- mm[, !(colnames(mm) %in% c('duration12'))]
    # mm <- mm[, !(colnames(mm) %in% c('duration6:snTRUE'))]
    # mm <- mm[, !(colnames(mm) %in% c('duration24:snTRUE'))]
    # mm <- mm[, !(colnames(mm) %in% c('duration2:snTRUE'))]
    if (T) {
      QR <- qr(mm)
      if (!ncol(mm) == QR$rank) {
        message('Not full rank')
        # qr.R(QR)
        # M <- strip_dimnames(round(qr.R(QR) != 0))
        M <- round(qr.R(QR) != 0)
        M <- matrix(as.logical(M), nrow = nrow(M))
        redundant_col <-
          colnames(mm)[which(upper.tri(M) & !M, arr.ind = T)[, 2]]
        # mm <- mm[, !(colnames(mm) %in% redundant_col)]
      }
    }
  } else {
    mm <- model.matrix(
      ~1 + sn_dilution + duration, data = sa
    )
  }

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_6493 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data

  mm <- model.matrix(~1 + stim_group * duration, data = sa)

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_6600 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data

  mm <- model.matrix(~1 + stim_group, data = sa)

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_6601 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)
  sa <- so@meta.data
  sa$stim_group <- relevel(sa$stim_group, 'Unexposed in vivo')

  if (!any(c('tnf_conc', 'ifn_conc') %in% colnames(sa))) {
    extra_cols <- tar_read(sc_6601_sample_annotation)[sa$condition_i, ] %>%
      dplyr::select(tnf_conc, ifn_conc)
    sa <- cbind(sa, extra_cols)
  }

  if (T) {
    sa$ifn_conc <- numeric2factor(sa$ifn_conc)
    sa$tnf_conc <- numeric2factor(sa$tnf_conc)
    sa$duration <- numeric2factor(sa$duration)

    mm <- model.matrix(~1 + tnf_conc * ifn_conc + duration + ifn_conc:duration,
      data = sa)
    mm <- mm[, apply(mm, 2, sum) > 0]
    mm <- mm[, !colnames(mm) %in% c('duration24', 'duration48')]
    mm <- mm[, !colnames(mm) %in% c('ifn_conc100:duration24')]
  } else {
    mm <- model.matrix(~1 + stim_group + duration, data = sa) %>%
      bind_cols(
        matrix(
          as.integer(with(sa, duration == 48 &
              stim_group == '100 ng/ml IFNy'))) %>%
        set_colnames('stim_group100 ng/ml IFNy:duration48')
      )
  }

  if (F) {
    # mm <- model.matrix(~1 + stim_group*duration, data = sa)
    mm <- model.matrix(~1 + stim_group + stim_group:duration, data = sa)
    # mm <- mm[, !stringr::str_detect(colnames(mm), 'Unexposed')]
    # mm <- mm[, !stringr::str_detect(colnames(mm), 'TNFa.*24')]
    mm <- mm[, colnames(mm) != 'stim_group10 ng/ml TNFa:duration24']
    mm <- mm[, colnames(mm) != 'stim_group10 ng/ml TNFa:duration48']
    mm <- mm[, colnames(mm) != 'stim_group100 ng/ml IFNy 10 ng/ml TNFa:duration48']
    mm <- mm[, colnames(mm) != 'stim_group100 ng/ml IFNy 10 ng/ml TNFa:duration24']
    mm <- mm[, colnames(mm) != 'stim_group100 ng/ml IFNy:duration24']
  }

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


#' 6601, due to it's uneven design, focus on stim group (sg)
#'
#'
limma_model_6601_sg <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)
  sa <- so@meta.data
  sa$stim_group <- relevel(sa$stim_group, 'Unexposed in vivo')

  mm <- model.matrix(~1 + stim_group * duration, data = sa)
  mm <- mm[, apply(mm, 2, sum) > 0]
  # mm <- mm[, !colnames(mm) %in% c('duration24', 'duration48')]
  mm <- mm[, !colnames(mm) %in% c('stim_group100 ng/ml IFNy:duration24')]

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_6434 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data
  stopifnot(rownames(sa) == colnames(d0))
  sa <- separate_duration(sa)
  # sa$sn_dilution <- stringr::str_replace(sa$stim_group, ' SN$', '')
  sa$ifn <- sa$ifn_duration > 0
  sa$tnf <- sa$tnf_duration > 0
  sa <- sa %>%
    dplyr::mutate(
      across(
        matches('duration|conc|dilution'),
        numeric2factor
    ))
  mm <- model.matrix(
    ~1 + tnf_conc + ifn_conc + duration + ifn:tnf:duration,
    data = sa
  )
  allowed_cols <- stringr::str_subset(colnames(mm), '^[^:]*$|:.*TRUE')
  mm <- mm[, allowed_cols]
  mm <- mm[, !(colnames(mm) %in% c('tnf_conc0.1', 'ifn_conc1'))]

  if (F) {
    QR <- qr(mm)
    ncol(mm) == QR$rank
    qr.R(QR)
    M <- strip_dimnames(round(qr.R(QR) != 0))
    M <- matrix(as.logical(M), nrow = nrow(M))
  }

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


limma_model_5029 <- function(so, down_scale = 1) {
  d0 <- limma_preproc(so = so, down_scale = down_scale)

  sa <- so@meta.data
  stopifnot(rownames(sa) == colnames(d0))
  if (T) {
    sa <- separate_duration(sa)
    # sa$sn_dilution <- stringr::str_replace(sa$stim_group, ' SN$', '')
    sa$ifn <- sa$ifn_duration > 0
    sa$tnf <- sa$tnf_duration > 0
    sa$sn <- sa$sn_dilution > 0
    sa <- sa %>%
      dplyr::mutate(
        across(
          matches('duration|conc|dilution'),
          numeric2factor
      ))
    mm <- model.matrix(
      ~1 + tnf_conc + ifn_conc + sn_dilution + duration +
        sn:duration + ifn:tnf:duration,
      data = sa
    )
    allowed_cols <- stringr::str_subset(colnames(mm), '^[^:]*$|:.*TRUE')
    mm <- mm[, allowed_cols]
    mm <- mm[, !(colnames(mm) %in%
      c('sn_dilution5e-05', 'tnf_conc0.1', 'ifn_conc1'))]

    if (F) {
      QR <- qr(mm)
      ncol(mm) == QR$rank
      qr.R(QR)
      M <- strip_dimnames(round(qr.R(QR) != 0))
      M <- matrix(as.logical(M), nrow = nrow(M))
    }
  } else {
  }

  efit <- estimate_limma(d0 = d0, mm = mm, span = .3)

  return(efit)
}


extract_terms <- function(
  efit,
  min_estimate = 1,
  max_p = .25,
  verbose = F,
  regex = '.*') {

  pacman::p_load('biobroom')

  t_dat <-
    tidy(efit) %>%
    dplyr::mutate(p.value.adj = p.adjust(p.value, 'fdr')) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(term = factor(term, levels = unique(term)))

  time_informative_genes <-
    t_dat %>%
    dplyr::filter(
      abs(estimate) >= min_estimate & p.value.adj <= max_p) %>%
    dplyr::filter(stringr::str_detect(term, regex)) %>%
    dplyr::select(term, gene, statistic, estimate, p.value.adj) %>%
    dplyr::group_by(term) %>%
    dplyr::arrange(p.value.adj, .by_group = TRUE) %>%
    dplyr::ungroup()

  if (verbose)
    print(table(time_informative_genes$term))

  return(time_informative_genes)
}


diagnose_mm_problems <- function(mm) {
  if (F) {
    colnames(mm)
    QR <- qr(mm)
    ncol(mm) == QR$rank
    strip_dimnames(round(qr.R(QR) != 0))
    # qr.Q(QR)
    # str(QR)
    # QR$qr
    # QR$pivot
    # QR$qraux
    # str(QR$qr)
    # mm
  }
  # abs(mm[, 'duration24'] - mm[, 'ifn_duration24']) %>%
  #   { .[. > 0] }
  # local({
  #   NR <- nrow(mm)
  #   attributes(mm) <- NULL
  #   matrix(mm, nrow = NR)
  # })
}


comparison_volcano <- function(
  experiment_i, experiment_j,
  res_i, res_j,
  fn_add = '', x_name = '', y_name = '') {

  if (maartenutils::null_dat(res_i) ||
      maartenutils::null_dat(res_j)) {
    return('')
  }
  # pacman::p_load('emojifont')

  # res_i$comp <- 'A'
  # res_j$comp <- 'B'
  # p_dat <- dplyr::bind_rows(res_i, res_j)
  p_dat <- dplyr::full_join(res_i, res_j, by = 'gene')
  p_dat <- p_dat %>%
    dplyr::mutate(max_estimate =
      pmax(estimate.x, estimate.y, na.rm = T)) %>%
    dplyr::mutate(min_p =
      pmin(p.value.adj.x, p.value.adj.y, na.rm = T)) %>%
    dplyr::mutate(estimate.x = ifelse(is.na(estimate.x), 0, estimate.x)) %>%
    dplyr::mutate(estimate.y = ifelse(is.na(estimate.y), 0, estimate.y))


  repl_NA <- function(x) {
    x[is.na(x)] <- FALSE 
    return(x)
  }
  stats <- 
    p_dat %>%
    dplyr::ungroup() %>%
    dplyr::filter(min_p <= .1) %>%
    dplyr::summarize(
      # X = sum(!is.na(term.x))/n(),
      # Y = sum(!is.na(term.y))/n()
      X = sum(repl_NA(p.value.adj.x <= .1))/n(),
      Y = sum(repl_NA(p.value.adj.y <= .1))/n()
    )

  mme <- quantile(p_dat$max_estimate, .99, na.rm = T)
  hl_points <- p_dat %>%
    dplyr::filter(max_estimate >= mme & min_p <= .1)

  p <- ggplot(p_dat, aes(x = estimate.x, y = estimate.y)) +
    ggrepel::geom_text_repel(
      min.segment.length = 0.01,
      show.legend = FALSE,
      size = 2, data = hl_points, mapping = aes(label = gene)
    ) +
    geom_point(alpha = .1, 
      mapping = aes(size = -log10(min_p))) +
    # geom_point(alpha = .1, shape = 'U+25C2') +
    # geom_point(alpha = .1, shape = 'U+25B8') +
    # geom_point(
    #   alpha = .1,
    #   mapping = aes(size = log10(1/p.value.adj.x)),
    #   # shape = 'U+25D6') +
    #   shape = '\u25D6') +
    # geom_point(
    #   alpha = .1,
    #   mapping = aes(size = log10(1/p.value.adj.y)),
    #   shape = '\u25D7') +
    # xlab(glue::glue('{x_name} - {100 *round(stats$X, 2)}%')) +
    # ylab(glue::glue('{y_name} - {100 *round(stats$Y, 2)}%')) +
    xlab(glue::glue('{x_name}')) +
    ylab(glue::glue('{y_name}')) +
    scale_size_continuous(name = '-log10(p)', range = c(0, 3))

  return(p)
}
