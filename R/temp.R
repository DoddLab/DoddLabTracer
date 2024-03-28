# load('~/Project/04_package/MetDNA2/data/lib_adduct_nl.rda')
# usethis::use_data(lib_adduct_nl)
#
# writexl::write_xlsx(lib_adduct_nl, path = '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/lib_adduct_nl_231003.xlsx', format_headers = FALSE)
#
# lib_adduct_nl_pos <- readxl::read_xlsx('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/lib_adduct_nl_modified_231003.xlsx', sheet = 1)
# lib_adduct_nl_neg <- readxl::read_xlsx('~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/lib_adduct_nl_modified_231003.xlsx', sheet = 2)
# lib_adduct_nl <- list('positive' = lib_adduct_nl_pos,
#                       'negative' = lib_adduct_nl_neg)
# usethis::use_data(lib_adduct_nl, overwrite = TRUE)
