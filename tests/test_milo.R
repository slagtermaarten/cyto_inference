coq_sc_so <- tar_read_raw('filtered_so_6369')
q_sc_so <- tar_read_raw('filtered_so_6489')

source('~/MirjamHoekstra/R/init.R')
sce <- perform_agg_neighbourhoods(
  so = q_sc_so,
  genelist = 'informativeV15',
  k = 100,
  d = 5
)


# sce <- tar_read(agg_neighbourhoods_6369_informativeV15_500_5)
sce <- tar_read(agg_neighbourhoods_6600_informativeV15_300_5)
sce <- tar_read(agg_neighbourhoods_6489_informativeV15_100_5)
if (is.null(miloR::nhoodDistances(sce))) {
  sce <- calcNhoodDistance(sce, d=5, reduced.dim = 'pca.corrected')
}

source('~/MirjamHoekstra/R/init.R')
DA_Nhoods <- test_milo_DA_Nhoods(sce, debug_f = T)

ref_so <- tar_read(bulk_comb_so)
NH_so <- tar_read(NH_so_6369_informativeV15_500_5_10)

# NH_so@assays
# GetAssayData(ref_so, 'counts', assay = 'SCT')[1:5, 1:5]
# GetAssayData(ref_so, 'counts', assay = 'SCT')[1:5, 1:5]
# GetAssayData(NH_so, 'counts', assay = 'RNA')[1:5, 1:5]
# SetAssayData(NH_so, 'counts', object = )
M <- as.matrix(GetAssayData(NH_so, 'counts', assay = 'RNA')) 
# NH_so@assays$SCT <- SetAssayData(NH_so, slot = 'counts',
#   new.data = M, assay = 'SCT')
# NH_so[['SCT']] <- SetAssayData(NH_so, slot = 'counts',
#   new.data = M, assay = 'SCT')
# NH_so <- SetAssayData(NH_so, slot = 'counts',
#   new.data = M, assay = 'SCT')
NH_so[["SCT"]] <- CreateAssayObject(counts = M)

source('~/MirjamHoekstra/R/init.R')
out <- extract_transact_Ms(
  reference_so = ref_so,
  query_so = NH_so,
  normalize_sample_list_params = list(
    genes = NULL,
    genelist = 'informativeV15',
    GDR_thresh = 0,
    norm_method = 'TMM',
    Z_scale = T
  )
)

Ms <- tar_read(NH_harmonized_Ms_comb_6600_informativeV15_30_3_10)
purrr::map(Ms, ~rowSums(.x))
