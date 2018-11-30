test_df <- tibble(tnf_conc = c(10, 10, 10), tnf_duration = c(2, 2, 12))
attr(test_df, 'name') <- 'test_name'
stopifnot(attributes(numerify_regressors(test_df))$name == 'test_name')

recover_garbled_string(
  query = '1.ng.ml.IFNy...6h',
  candidates = c('1 ng/ml IFNy 6h', '1 ng/ml IFNy 2h'),
  N_return = 2L
)

find_var_part()
find_var_part('aaa')
