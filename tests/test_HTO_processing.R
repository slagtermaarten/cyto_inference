## These are the first few cells of exp 6600
M <- structure(c(0.0773943949008814, 0.0394457472014829, 0.216243295589283,
2.47513520618877, 0.13174661495026, 1.84227768712677, 2.79384084968956,
0.0352220592398676, 1.91276725901465, 0.102149564871066, 3.69702430647829,
2.59822395342984, 0.0534918590650635, 2.74831366028125, 0.198681675471593,
0.104267044275215, 0.325439328926619, 0.198681675471593, 0.136051873936147,
0.547373516670556, 0.362713460293686, 0.547373516670556, 0.459299863200042,
0.362713460293686, 0.335502364203314, 0, 0.0950016860565839,
0.139319291221398, 0.437050270574602, 2.51859535922403, 0, 0.0969743037303906,
2.71156563598919, 0.0969743037303906, 0.0969743037303906, 0.476850362843634,
3.08146431506563, 0.330935195454842, 0.282838840408089, 3.01958168578059,
3.82250806720944, 2.94677309213488, 0.394118325449049, 0.676078357421634,
0, 0.394118325449049, 0.394118325449049, 0.394118325449049), .Dim = c(6L,
8L), .Dimnames = list(c("6600_AAACCCAAGTGTTCCA", "6600_AAACCCAAGTTGTCGT",
"6600_AAACCCACAGTAGTGG", "6600_AAACCCAGTTTGGAGG", "6600_AAACCCATCACTACGA",
"6600_AAACCCATCGATACTG"), c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5",
"HTO6", "HTO7", "HTO8")))

# dput(hash_tag_annotation)
hash_tag_annotation <- 
  structure(list(HT1 = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
      HT2 = c(4, 5, 6, 7, 8, 4, 5, 6, 8, 5, 6, 7, 8), 
      stim_group = c("10 ng/ml TNFa",
      "10 ng/ml TNFa", "10 ng/ml TNFa", "100 ng/ml IFNy 10 ng/ml TNFa",
      "100 ng/ml IFNy 10 ng/ml TNFa", "Unstimulated in vitro",
      "100 ng/ml IFNy", "Unstimulated in vitro", "Unstimulated in vitro",
      "Unstimulated in vitro", "Unstimulated in vitro", "100 ng/ml IFNy",
      "100 ng/ml IFNy"), duration = structure(c(2L, 2L, 2L, 2L,
      2L, 4L, 5L, 5L, 4L, 2L, 2L, 4L, 4L), .Label = c("2", "6",
      "12", "24", "48", "Unknown"), class = "factor"), mouse = c("6601_1",
      "6601_2", "6601_3", "6601_4", "6601_5", "6601_11", "6601_12",
      "6601_13", "6601_10", "6601_6", "6601_7", "6601_8", "6601_9"
      ), condition_i = 1:13), row.names = c(NA, -13L), 
  class = c("data.table", "data.frame"))

sc_HT_diagnostics <- 
  find_matching_double_barcodes(
    M = M, 
    hash_tag_annotation = hash_tag_annotation, 
    debug = F,
    ncores = 1)

testthat::expect_equal(
  sc_HT_diagnostics[, c('HT1', 'HT2')],
  structure(list(HT1 = c(2, 3, 2, 1, 2, 2), HT2 = c(7, 8, 6, 7,
  7, 7)), class = "data.frame", row.names = c(NA, -6L)))


# saveRDS(M, file.path('6600_proc_HTO.rds'))
M <- readRDS(file.path('6600_proc_HTO.rds'))
# dput(hash_tag_annotation)
assign('sa', tar_read(sc_6600_sample_annotation), envir = .GlobalEnv)
hash_tag_annotation <- sa

source('R/HTO_processing.R')
tmp <- 
  find_matching_double_barcodes(
    M = M[1:5, , drop = F], 
    hash_tag_annotation = hash_tag_annotation, 
    debug = F)

