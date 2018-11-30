list_genesets <- function() {
  paths <- 
    c(file.path(data_raw_dir), file.path(data_raw_dir, 'gene_lists'))
  purrr::map_dfr(paths, 
    maartenutils::gen_file_overview, pat = '.*genes.txt', 
    include_full = T) %>%
    dplyr::arrange(desc(mtime)) %>%
    dplyr::mutate(N_genes = map_dbl(full_fn, 
        ~length(read_geneset_fn(.x)))) %>%
    dplyr::mutate(name = basename(full_fn)) %>%
    # dplyr::select(-full_fn)
    { . }
}


read_geneset_fn <- function(fn) {
  genes <- readLines(con = fn)
  genes <- genes[!grepl('^#', genes)]
  genes
}


read_geneset <- function(name) {
  if (is.null(name) || is.na(name) || 
      (is.character(name) && name == 'NA')) {
    return(NULL)
  }
  # list.files(data_raw_dir)
  potential_fns <- 
    c(file.path(data_raw_dir, 
        glue::glue('{name}.txt')),
      file.path(data_raw_dir, 'gene_lists', 
        glue::glue('{name}.txt')),
      file.path(data_raw_dir, 
        glue::glue('{name}_genes.txt')),
      file.path(data_raw_dir, 'gene_lists', 
        glue::glue('{name}_genes.txt')))
  # find_character_diff(
  #   potential_fns[4],
  #   flatten_chr(list_genesets()[1, 'full_fn'])
  # )
  fn <- potential_fns[file.exists(potential_fns)]
  if (length(fn) > 1) {
    rlang::warn(paste('Found:', paste(fn, collapse = ', ')))
    fn <- fn[1]
  }
  if (length(fn) == 0) stop('geneset not found')
  genes <- read_geneset_fn(fn)
  return(genes)
}
# read_geneset('informativeV15')


mirjam_genes_220516 <- list(
  TNFa_mono = "Cd34
    Map4k4
    Pdgfb
    Ets2
    Smad7
    Orai2
    Abhd17c
    Tnnt2
    Sema6d
    Dclk1
    Spsb1
    Col7a1
    Mmp12
    Nfkbia
    Kcnk10
    Slc7a8
    Nfkbie
    Syt7
    Eif6
    Vcam1
    Tlr2
    Cxcl1
    Clec4e
    Relb
    Birc3
    Tagln2
    Stk10
    Itgb3 
    Hid1
    Trp53
    Odc1
    Itga2
    A4galt
    Cd82
    Ass1
    Pnrc1
    Anpep
    Mmp1b
    Mmp10
    Maged2
    Il13ra2
    Ern1
    Sh3pxd2b
    Zhx2
    Ifnar1
    5730559C18Rik
    Ackr3 
    Map2k6
    Rin2
    Micall2
    Xylt1
    Cdh5
    Adamts7
    Nbeal2
    Tnip1
    Traf1
    Ptges
    Rnd1
    Ifngr2
    Plpp3
    Mmp3
    Psmd10
    Ifngr1
    Ubtd2
    Ppic
    Tmem176b
    Smoc2
    Gjb4
    Efnb1
    Zeb1
    Sh3pxd2a
    Notch2
    Gpr149
    Sh2d5
    Penk
    Apbb2
    Acsbg1
    Abcc4 
    Denr
    Gprc5a
    " %>% line_split(),

  TNFa_mono_doubt = "Rgs16
    Neurl3
    Rnf149
    Angptl4
    Rela
    Nfkb2
    Tnip3
    Plek
    Ccbe1
    Slc43a3
    Sema4b
    Zswim4
    Acpp
    Timp1
    Gm9115
    Sele
    Rap2a
    Npc1
    Arap3
    Tnfaip3
    Stx11
    Plau
    Prkcd
    Glis3
    Bcl10
    Cemip
    Rab11fip1
    Rffl
    Il6st
    Plaur
    Nrg1
    Pdgfa
    Lrrc8b
    Fam84a
    Mfhas1" %>% line_split(),

  IFNy_mono = 'Tubb3
    Fgl2
    Mtfr2
    Gbp11
    C4b
    Kctd3
    Nfyb
    2610507B11Rik
    Atxn7l1
    Fst
    Crebrf
    Fam171b
    Slit2
    Tnfrsf1a
    Clic5
    Bambi
    Fhod3
    Il33
    Batf2
    Gbp10
    Ppp1r3b
    Ppp1r26
    Tigd2
    Gm12185
    Twist1
    Pfkp
    Ly6c1
    Tmem173 
    Phactr4
    Loxl3
    Cnep1r1
    Kitl
    Fam107b
    Apol6
    Socs1
    Crem
    Gm4841
    Prune2
    Lmo4
    Acsl1
    Fos
    # Tinagl1
    Rock1
    Lrrtm2
    Gclm
    Aebp2
    Cx3cl1
    Irf8
    Pcdh18
    Irak2
    Sema3e
    Styk1
    Pdzrn3
    Mpp1' %>% line_split(),

  IFNy_mono_doubt = '
    Nfil3
    Amigo2
    Greb1l
    Sat1
    Gbp4
    Rgs12
    Asap3
    Nudt9
    C2cd5
    Chic1
    Rab3ip
    Zfp367
    Slc4a11
    Casp8ap2
    Htt
    Otud1
    Cldn12
    Socs3
    Wars 
    Trappc11
    Casp12
    ' %>% line_split(),

  IFNy_mono_synergy = '
    Nos2
    Slamf8
    Rhob2
    Serpina3g
    Serpina3f
    Ubd
    Cxcl9
    Gbp8
    C1s1
    ' %>% line_split()
)

mirjam_genes_220518 <- list(
  TNFa_early = "
    Map4k4
    Pdgfb
    Ets2
    Smad7
    Orai2
    Nfkbia
    Nfkbie
    Vcam1
    Cxcl1
    Zhx2
    Ifngr1
    Ubtd2
    Gjb4
  " %>% line_split(),

  TNFa_late = "
      Tnnt2
      Sema6d
      Dclk1
      Spsb1
      Odc1
      Itga2
      A4galt
      Ass1
      Cd82
      Mmp1b
      Mmp10
      Maged2
      Il13ra2
      Map2k6
      Rin2
      Xylt1
      Cdh5
      Adamts7
      Nbeal2
      Gpr149
  " %>% line_split(),

  synergy_plateau = "
    Cxcl9
    Ubd
  " %>% line_split(),

  synergy_late = "
    Serpina3f
    Gbp8
    C1s1
  " %>% line_split(),

  IFNy_early = "
    Fgl2
    Atxn7l1
    Fst
    Tnfsrf1a
    Kitl
    Socs1
    Crem
    Lmo4
    Rock1
    Lrrm2
    Gclm
    Aebp2
    Cx3xl1
    Irf8
    Irak2
    Pdzrn3" %>% line_split(),

  IFNy_late = "
      Gbp11
      C4b
      Fhod3
      Il33
      Batf2
      Gbp10
      Ppp1r3b
      Ppp1r26
      Tigd2
      Gm12185
      Twist1
      Ly6c1
      Phactr4
      Loxl3
      Cnep1r1
      Tmem200a
  " %>% line_split()
)

mirjam_genes_220519 <- list(
  synergy_less_strict = "
    Cxcl9
    Ubd
    Serpina3f
    Serpina3g
    Gbp8
    C1s1
    Nos2
    Slamf8
    Gbp4 
    Asap3 
    Phyh 
    Pla2r1 
    Gem 
    Cdh17 
    Arhgdib 
    Trim30b 
    Il15 
    Ifi209 
    Serpinb2 
    Ccl7 
    Ccl2 
    Ptch1 
    Gsdmd 
    Pim1 
    H2−Q7 
    H2−Q6 
    Fyn 
    Prl2c3 
    Htr2a
    Icam1
    Tmg2
    Slc39a1
    Fas
    Rnd1
    Phlda
    Enc1
    Il15ra
    Ccrl2
    Casp4
  " %>% line_split()
)

mouse_genesets <- 
  sapply(ls(pattern = 'mirjam_genes_2205'), get) %>%
  purrr::flatten() %>%
  { .[!stringr::str_detect(names(.), 'doubt')] }

if (F) {
  mouse_genesets  %>%
    { .[!stringr::str_detect(names(.), 'mono_synergy')] } %>%
    { .[stringr::str_detect(names(.), 'mono|less_strict')] } %>%
    # map_dbl(length) %>%
    # sum()
    unlist() %>%
    unique() %>%
    length()
}
