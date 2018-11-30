plotting_targets <- list(
  ### ------------   Plotting related   ------------ ###
  tar_target(stim_cols, {
    old_colors <- gen_cyto_inf_col_scale()
    cs <-
      dplyr::bind_rows(
        sample_annotation_exp6434,
        sample_annotation_exp5029
      ) %>%
      gen_cyto_inf_col_scale_from_obs()
    old_colors[names(cs)] <- cs
    return(old_colors)
 })
)
