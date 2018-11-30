extract_gene_stats <- function(efit, M_cpm, LL = .5, UL = .2, SD = .8,
  FD = 4, p_u_thresh = 1e-1, p_l_thresh = 1e-1) {
  all_terms_l <- 
    extract_terms(efit,
      max_p = 1,
      min_estimate = 0
    )

  gene_max_e <- apply(M_cpm, 1, max)

  # all_terms_l %>%
  #   dplyr::filter(gene %in% 'Traf1')

  out <- 
    all_terms_l %>%
    dplyr::group_by(gene) %>%
    # dplyr::filter(p.value.adj[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'] <= .25) %>%
    # dplyr::filter(estimate[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'] > 0) %>%
    # dplyr::filter(estimate[term == 'stim_group10 ng/ml TNFa'] > 0) %>%
    # dplyr::filter(estimate[term == 'stim_group100 ng/ml IFNy'] > 0) %>%
    dplyr::summarize(
      tnfa = statistic[term == 'stim_group10 ng/ml TNFa'],
      tnfa_p = p.value.adj[term == 'stim_group10 ng/ml TNFa'],
      ifny = statistic[term == 'stim_group100 ng/ml IFNy'],
      ifny_p = p.value.adj[term == 'stim_group100 ng/ml IFNy'],
      tgfb = statistic[term == 'stim_group10 ng/ml TGFb'],
      tgfb_p = p.value.adj[term == 'stim_group10 ng/ml TGFb'],
      ifna = statistic[term == 'stim_group5000 U/ml IFNA1'],
      ifna_p = p.value.adj[term == 'stim_group5000 U/ml IFNA1'],
      combo = statistic[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'],
      combo_p = p.value.adj[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'],
      synergy = 
        (statistic[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'] -
          statistic[term == 'stim_group100 ng/ml IFNy'] - 
          statistic[term == 'stim_group10 ng/ml TNFa']) /
        (statistic[term == 'stim_group100 ng/ml IFNy'] + 
          statistic[term == 'stim_group10 ng/ml TNFa'])
    ) %>%
    dplyr::mutate(synergy = ifelse(combo_p <= .1 & combo > 0, synergy, 0)) %>%
    # dplyr::summarize(
    #   synergy_t = 
    #     (t[term == 'stim_group100 ng/ml IFNy 10 ng/ml TNFa'] -
    #       estimate[term == 'stim_group100 ng/ml IFNy'] - 
    #       estimate[term == 'stim_group10 ng/ml TNFa']) /
    #     (estimate[term == 'stim_group100 ng/ml IFNy'] + 
    #      estimate[term == 'stim_group10 ng/ml TNFa'])
    # ) %>%
    # dplyr::arrange(desc(synergy)) %>%
    dplyr::left_join(tibble::enframe(gene_max_e, 'gene', 'max_e')) %>%
    dplyr::mutate(min_p = apply(across(matches('_p$')), 1, min)) %>%
    dplyr::mutate(max_statistic = apply(across(c(tnfa, ifny, combo)), 1, max)) %>%
    dplyr::summarize(across(), across(c(tnfa, ifny, tgfb), 
      list(
        'neg' = ~abs(.x) <= UL, 
        'pos' = ~.x >= LL, 
        ## e2c := estimate 2 combo, i.e., does the mono
        ## response differ from the combo response?
        'e2c' = ~abs(.x - combo) <= .1)
      )
    ) %>%
    ## _l stands for 'loose', i.e. less stringent criteria
    dplyr::summarize(across(), across(c(tnfa, ifny, tgfb), 
        list(
          'neg_l' = ~.x <= UL, 
          'pos_l' = ~.x >= LL, 
          ## e2c := estimate 2 combo, i.e., does the mono
          ## response differ from the combo response?
          'e2c_l' = ~abs(.x - combo) <= .2)
        )) %>%
    dplyr::select(-tgfb_e2c, -tgfb_e2c_l) %>%
    dplyr::mutate(ifny_mr = ifny_pos & ifny_e2c & tnfa_neg & ifny_neg) %>%
    dplyr::mutate(tnfa_mr = tnfa_pos & tnfa_e2c & ifny_neg & tgfb_neg) %>%
    dplyr::mutate(tgfb_mr = tgfb_pos & tnfa_neg & ifny_neg) %>%
    # dplyr::mutate(tnfa_mr_l = 
    #   (tnfa_pos_l & (ifny_neg_l | (tnfa - ifny > SD) | tnfa / ifny >= FD)) & 
    #   TRUE
    #   # & (tnfa_e2c_l | synergy < 0)
    # ) %>%
    dplyr::mutate(tnfa_mr_l = 
      (
        (
          (tnfa_pos_l | (tnfa_p <= p_u_thresh)) & 
          (ifny_neg_l | (ifny_p > p_l_thresh))
        ) |
        (tnfa - ifny > SD) | 
        (tnfa / ifny >= FD)
      ) &
      TRUE
      # !tnfa_pos_l & (ifny_e2c_l | synergy < 0)
    ) %>%
    dplyr::mutate(ifny_mr_l = 
      (
        (
          (ifny_pos_l | (ifny_p <= p_u_thresh)) & 
          (tnfa_neg_l | (tnfa_p > p_l_thresh))
        ) |
        (ifny - tnfa > SD) | 
        (ifny / tnfa >= FD)
      ) &
      TRUE
      # !tnfa_pos_l & (ifny_e2c_l | synergy < 0)
    ) %>%
    dplyr::mutate(tgfb_mr_l = tgfb_pos_l & tnfa_neg_l & ifny_neg_l) %>%
    dplyr::arrange(desc(max_statistic)) %>%
    { . }

  return(out)
}
