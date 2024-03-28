# #
#
# library(tidyverse)
# library(DoddLabTracer)
#
# find_intemidates(peak_table_unlabel = '230830_Csp_YgeX_WT_UA.xlsx',
#                  peak_table_label = '230830_Csp_YgeX_WT_13CUA.xlsx',
#                  path = '~/Project/00_Uric_Acid_project/Data/20230912_example/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 20)
#
#
#
#
# # peak_table <- 'ssnA_UA_48h_area.txt'
# # peak_table_label <- 'ssnA_13CUA_48h_area.txt'
# # path = '~/Project/00_Uric_Acid_project/Data/230928_isotope_tracing/ssnA/'
# # control_group = "WT_UA"
# # case_group = "ssnA_UA"
# #
# find_intemidates(peak_table_unlabel = 'ssnA_UA_48h_area.txt',
#                  peak_table_label = 'ssnA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess/ssnA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ssnA_UA', 'ssnA_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE)
#
#
# find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
#                  peak_table_label = 'hyuA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess/hyuA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('hyuA_UA', 'hyuA_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE)
#
#
# # find_intemidates(peak_table_unlabel = 'xdh_UA_48h_area.txt',
# #                  peak_table_label = 'xdh_13CUA_48h_area.txt',
# #                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess/xdh/',
# #                  control_group = c("WT_UA", "WT_13CUA"),
# #                  case_group = c('xdh_UA', 'xdh_13CUA'),
# #                  mz_tol = 10,
# #                  rt_tol = 0.05,
# #                  p_value_cutoff = 0.05,
# #                  fold_change_cutoff = 10,
# #                  p_adjust = FALSE)
#
#
# find_intemidates(peak_table_unlabel = 'ygeW_UA_48h_area.txt',
#                  peak_table_label = 'ygeW_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess/ygeW/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeW_UA', 'ygeW_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE)
#
#
# find_intemidates(peak_table_unlabel = 'ygeY_UA_48h_area.txt',
#                  peak_table_label = 'ygeY_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess/ygeY/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeY_UA', 'ygeY_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE)
#
#
# find_intemidates(peak_table_unlabel = 'ygfK_UA_48h_area.txt',
#                  peak_table_label = 'ygfK_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/230928_isotope_tracing/ygfK/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygfK_UA', 'ygfK_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10)
#
#
#
# ################################################################################
#
# # ssnA
# find_intemidates(peak_table_unlabel = 'ssnA_UA_48h_area.txt',
#                  peak_table_label = 'ssnA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ssnA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ssnA_UA', 'ssnA_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygfK
# find_intemidates(peak_table_unlabel = 'ygfK_UA_48h_area.txt',
#                  peak_table_label = 'ygfK_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygfK/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygfK_UA', 'ygfK_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # HyuA
# find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
#                  peak_table_label = 'hyuA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/hyuA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('hyuA_UA', 'hyuA_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygeW
# find_intemidates(peak_table_unlabel = 'ygeW_UA_48h_area.txt',
#                  peak_table_label = 'ygeW_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeW/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeW_UA', 'ygeW_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygeY
# find_intemidates(peak_table_unlabel = 'ygeY_UA_48h_area.txt',
#                  peak_table_label = 'ygeY_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeY_UA', 'ygeY_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
# # load('~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/00_tracer_result/00_intermediate_data/list_peak_group_annotation.RData')
# # load('~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/00_tracer_result/00_intermediate_data/list_peak_group_annotation_concised.RData')
# # load('~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/00_tracer_result/00_intermediate_data/list_peak_group_annotation_merge.RData')
# #
# # list_peak_group_annotation$`131.0452@6.438`@peak_list_annotated
#
#
# # ygeX
# find_intemidates(peak_table_unlabel = 'ygeX_UA_48h_area.txt',
#                  peak_table_label = 'ygeX_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeX/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)

# ################################################################################
# # # ygeY
# load('~/Project/00_Uric_Acid_project/Data/20231108_isotope_tracing_analysis/ygeY/00_tracer_result/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('~/Project/00_Uric_Acid_project/Data/20231108_isotope_tracing_analysis/ygeY/00_tracer_result/00_intermediate_data/list_peak_group_annotation_merge.RData')
#
# temp <- list_peak_group_annotation_concised$`148.0715@6.437`
# temp <- list_peak_group_annotation_merge$`148.0715@6.437`
# plotPseudoMs1Spec(temp)
#
# temp <- list_peak_group_annotation_merge$`190.082@4.88`
# plotPseudoMs1Spec(temp)


################################################################################

# ################################################################################
#
# # ssnA
# find_intemidates(peak_table_unlabel = 'ssnA_UA_48h_area.txt',
#                  peak_table_label = 'ssnA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ssnA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ssnA_UA', 'ssnA_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = FALSE)

#
# # ygfK
# find_intemidates(peak_table_unlabel = 'ygfK_UA_48h_area.txt',
#                  peak_table_label = 'ygfK_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygfK/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygfK_UA', 'ygfK_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = FALSE)
#
#
# # HyuA
# find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
#                  peak_table_label = 'hyuA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/hyuA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('hyuA_UA', 'hyuA_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygeW
# find_intemidates(peak_table_unlabel = 'ygeW_UA_48h_area.txt',
#                  peak_table_label = 'ygeW_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeW/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeW_UA', 'ygeW_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygeY
# find_intemidates(peak_table_unlabel = 'ygeY_UA_48h_area.txt',
#                  peak_table_label = 'ygeY_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeY/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeY_UA', 'ygeY_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#
# # ygeX
# find_intemidates(peak_table_unlabel = 'ygeX_UA_48h_area.txt',
#                  peak_table_label = 'ygeX_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20231002_isotope_tracing_reprocess_copy/ygeX/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
#
#


# ################################################################################
#
# raw_data_unlabel <- modify_msdial_table(peak_table = 'ygeY_UA.txt',
#                                         path = '~/Project/00_Uric_Acid_project/Data/20231115_isotope_tracing_mutiple_timepoints',
#                                         control_group = 'WT_UA',
#                                         case_group = 'ygeY_UA')

################################################################################
#
# find_intemidates(peak_table_unlabel = 'hyuA_UA_48h_area.txt',
#                  peak_table_label = 'hyuA_13CUA_48h_area.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20240319_isotope_tracing_analysis/hyuA/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('hyuA_UA', 'hyuA_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)


################################################################################
#
#
# find_intemidates(peak_table_unlabel = '231227_ygeW_UA.txt',
#                  peak_table_label = '231227_ygeW_13CUA.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20240319_isotope_tracing_analysis/ygeW/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeW_UA', 'ygeW_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)

# # ygeX
# find_intemidates(peak_table_unlabel = '231227_ygeX_UA.txt',
#                  peak_table_label = '231227_ygeX_13CUA.txt',
#                  path = '~/Project/00_Uric_Acid_project/Data/20240319_isotope_tracing_analysis/ygeX/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#                  polarity = 'positive',
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 10,
#                  p_adjust = FALSE,
#                  is_recognize_adducts = TRUE)
