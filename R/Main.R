################################################################################
#' @title find_intemidates
#' @author Zhiwei Zhou
#' @param peak_table_unlabel unlabeled peak table
#' @param peak_table_label labeled peak table
#' @param path working directory. Default: "."
#' @param control_group label of control groups. e.g. c("WT_UA", "WT_13CUA")
#' @param case_group label of case groups. e.g. c('ygeX_UA', 'ygeX_13CUA')
#' @param mz_tol Default: 10 ppm
#' @param rt_tol Default: 0.05 min
#' @param p_value_cutoff Default: 0.05
#' @param fold_change_cutoff Default: 20
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' find_intemidates(peak_table_unlabel = '230830_Csp_YgeX_WT_UA.xlsx',
#'                  peak_table_label = '230830_Csp_YgeX_WT_13CUA.xlsx',
#'                  path = '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/',
#'                  control_group = c("WT_UA", "WT_13CUA"),
#'                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#'                  mz_tol = 10,
#'                  rt_tol = 0.05,
#'                  p_value_cutoff = 0.05,
#'                  fold_change_cutoff = 20)
#' }

# find_intemidates(peak_table_unlabel = '230830_Csp_YgeX_WT_UA.xlsx',
#                  peak_table_label = '230830_Csp_YgeX_WT_13CUA.xlsx',
#                  path = '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/',
#                  control_group = c("WT_UA", "WT_13CUA"),
#                  case_group = c('ygeX_UA', 'ygeX_13CUA'),
#                  mz_tol = 10,
#                  rt_tol = 0.05,
#                  p_value_cutoff = 0.05,
#                  fold_change_cutoff = 20)

# path <- '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/'
# peak_table_unlabel <- '230830_Csp_YgeX_WT_UA.xlsx'
# peak_table_label <- '230830_Csp_YgeX_WT_13CUA.xlsx'
# control_group = c("WT_UA", "WT_13CUA")
# case_group = c('ygeX_UA', 'ygeX_13CUA')
# mz_tol = 10
# rt_tol = 0.05
# p_value_cutoff = 0.05
# fold_change_cutoff = 20
#
# path <- '~/Project/00_Uric_Acid_project/Data/230928_isotope_tracing/ssnA/'
# peak_table_unlabel <- 'ssnA_UA_48h_area.txt'
# peak_table_label <- 'ssnA_13CUA_48h_area.txt'
# control_group = c("WT_UA", "WT_13CUA")
# case_group = c('ssnA_UA', 'ssnA_13CUA')
# mz_tol = 10
# rt_tol = 0.05
# p_value_cutoff = 0.05
# fold_change_cutoff = 20

# path <- '~/Project/00_Uric_Acid_project/Data/230928_isotope_tracing/ygeW/'
# peak_table_unlabel <- 'ygeW_UA_48h_area.txt'
# peak_table_label <- 'ygeW_13CUA_48h_area.txt'
# control_group = c("WT_UA", "WT_13CUA")
# case_group = c('ygeW_UA', 'ygeW_13CUA')
# mz_tol = 10
# rt_tol = 0.05
# p_value_cutoff = 0.05
# fold_change_cutoff = 10


# path <- '~/Project/00_Uric_Acid_project/Data/20231003_isotope_tracing_development_hyuA/hyuA/'
# peak_table_unlabel <- 'hyuA_UA_48h_area.txt'
# peak_table_label <- 'hyuA_13CUA_48h_area.txt'
# polarity <- 'positive'
# control_group = c("WT_UA", "WT_13CUA")
# case_group = c('hyuA_UA', 'hyuA_13CUA')
# mz_tol = 10
# rt_tol = 0.05
# p_value_cutoff = 0.05
# fold_change_cutoff = 10
# p_adjust = FALSE


find_intemidates <- function(peak_table_unlabel,
                             peak_table_label,
                             sample_info = NULL,
                             path = '.',
                             polarity = c('positive', 'negative'),
                             control_group = c("WT_UA", "WT_13CUA"),
                             case_group = c('ygeX_UA', 'ygeX_13CUA'),
                             mz_tol = 10,
                             rt_tol = 0.05,
                             p_value_cutoff = 0.05,
                             p_adjust = TRUE,
                             fold_change_cutoff = 20,
                             is_recognize_adducts = TRUE) {
  # browser()

  polarity <- match.arg(polarity)

  message(crayon::blue('1. Reading the peak tables... \n'))
  raw_data_unlabel <- modify_msdial_table(peak_table = peak_table_unlabel,
                                          path = path,
                                          control_group = control_group[1],
                                          case_group = case_group[1])


  raw_data_label <- modify_msdial_table(peak_table = peak_table_label,
                                        path = path,
                                        control_group = control_group[2],
                                        case_group = case_group[2])

  if (p_adjust) {
    raw_data_unlabel <- raw_data_unlabel %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        q_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(q_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change))

    raw_data_label <- raw_data_label %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        q_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(q_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change))
  } else {
    raw_data_unlabel <- raw_data_unlabel %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        p_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(p_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change))

    raw_data_label <- raw_data_label %>%
      dplyr::mutate(increase_label = dplyr::case_when(
        p_values <= p_value_cutoff & fold_change >= fold_change_cutoff ~ 'signif_increase',
        !(p_values < p_value_cutoff & fold_change > fold_change_cutoff) ~ 'unsignif')) %>%
      dplyr::mutate(log2_fc = log2(fold_change))
  }

  dir.create(file.path(path, '00_tracer_result', '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
  save(raw_data_unlabel, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'raw_data_unlabel.RData'))
  save(raw_data_label, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'raw_data_label.RData'))

# vocano plot
if (p_adjust){
  temp_plot1 <- ggplot2::ggplot(raw_data_unlabel) +
    ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(q_values),
                                     colour = increase_label)) +
    ggrepel::geom_text_repel(
      data = subset(raw_data_unlabel, increase_label == 'signif_increase'),
      ggplot2::aes(x = log2_fc, y = -log10(q_values), label = id),
      # position = position_jitter(seed = 1),
      box.padding = 1
    ) +
    ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
    ggplot2::ggtitle('Unlabeled UA - C. sporogenes ') +
    ZZWTheme()

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot1,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_unlabeled.pdf'),
                  width = 8, height = 6)


  temp_plot2 <- ggplot2::ggplot(raw_data_label) +
    ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(q_values), colour = increase_label)) +
    ggrepel::geom_text_repel(
      data = subset(raw_data_label, increase_label == 'signif_increase'),
      ggplot2::aes(x = log2_fc, y = -log10(q_values), label = id),
      # position = position_jitter(seed = 1),
      box.padding = 1
    ) +
    ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
    ggplot2::ggtitle('13C-Labeled UA - C. sporogenes') +
    ZZWTheme()

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot2,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_labeled.pdf'),
                  width = 8, height = 6)
} else {
  temp_plot1 <- ggplot2::ggplot(raw_data_unlabel) +
    ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(p_values),
                                     colour = increase_label)) +
    ggrepel::geom_text_repel(
      data = subset(raw_data_unlabel, increase_label == 'signif_increase'),
      ggplot2::aes(x = log2_fc, y = -log10(p_values), label = id),
      # position = position_jitter(seed = 1),
      box.padding = 1
    ) +
    ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
    ggplot2::ggtitle('Unlabeled UA - C. sporogenes ') +
    ZZWTheme()

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot1,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_unlabeled.pdf'),
                  width = 8, height = 6)


  temp_plot2 <- ggplot2::ggplot(raw_data_label) +
    ggplot2::geom_point(ggplot2::aes(x = log2_fc, y = -log10(p_values), colour = increase_label)) +
    ggrepel::geom_text_repel(
      data = subset(raw_data_label, increase_label == 'signif_increase'),
      ggplot2::aes(x = log2_fc, y = -log10(p_values), label = id),
      # position = position_jitter(seed = 1),
      box.padding = 1
    ) +
    ggplot2::geom_vline(xintercept = log2(fold_change_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed") +
    ggplot2::ggtitle('13C-Labeled UA - C. sporogenes') +
    ZZWTheme()

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot2,
                  filename = file.path(path, '00_tracer_result', 'vocano_plot_labeled.pdf'),
                  width = 8, height = 6)
}


  if (p_adjust) {
    feature_sig_unlabel <- raw_data_unlabel %>%
      dplyr::filter(q_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
    feature_sig_label <- raw_data_label %>%
      dplyr::filter(q_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
  } else {
    feature_sig_unlabel <- raw_data_unlabel %>%
      dplyr::filter(p_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
    feature_sig_label <- raw_data_label %>%
      dplyr::filter(p_values <= p_value_cutoff & fold_change >= fold_change_cutoff)
  }

  message(crayon::yellow('Unlabeled group: There are ', nrow(feature_sig_unlabel),
                         'features enriched ( p_value <=', p_value_cutoff,
                         '& fold-change >=', fold_change_cutoff, ')\n'))
  message(crayon::yellow('Labeled group: There are ', nrow(feature_sig_label),
                         'features enriched ( p_value <=', p_value_cutoff,
                         '& fold-change >=', fold_change_cutoff, ')\n'))

  if (is_recognize_adducts) {
    cat('\n')
    message(crayon::blue('1.1. Recognize adducts, neutral losses, in-source fragements of enriched peaks ... \n'))
    peak_table <- raw_data_unlabel %>%
      dplyr::rename('name' = id) %>%
      dplyr::select(name, mz, rt, dplyr::starts_with(control_group[1]), dplyr::starts_with(case_group[1]))

    peak_table_signif <- raw_data_unlabel %>%
      dplyr::filter(increase_label == 'signif_increase') %>%
      dplyr::rename('name' = id) %>%
      dplyr::select(name, mz, rt, dplyr::starts_with(control_group[1]), dplyr::starts_with(case_group[1]))


    dir.create(file.path(path, '00_tracer_result', '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
    save(peak_table, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'peak_table.RData'))
    save(peak_table_signif, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'peak_table_signif.RData'))

    list_peak_group_annotation_merge <- recognize_rela_peak(peak_table = peak_table,
                                                            peak_table_signif = peak_table_signif,
                                                            path_dir = path,
                                                            polarity = polarity,
                                                            raw_data_folder = case_group[1],
                                                            tol_mz = mz_tol,
                                                            tol_rt = 0.05,
                                                            cutoff_ssc = 0.8,
                                                            cutoff_ppc = 0.8)
    rm(peak_table, peak_table_signif);gc()
    result_peak_recognization_unlabel <- generate_recoginize_table(list_peak_group_annotation_merge = list_peak_group_annotation_merge)

    temp_num1 <- nrow(feature_sig_unlabel)
    temp_num2 <- length(unique(result_peak_recognization_unlabel$base_peak))

    message(crayon::green(temp_num1 - temp_num2, 'features are merged through recognizing adducts/neutral_loss/in_source_fragments\n'))
    message(crayon::green(temp_num2, 'unlabeled features used to extract pairs','\n'))
    rm(temp_num1, temp_num2);gc()

    feature_sig_unlabel <- feature_sig_unlabel %>% dplyr::filter(id %in% result_peak_recognization_unlabel$base_peak)

    save(feature_sig_unlabel, file = file.path(path, '00_tracer_result', '00_intermediate_data', 'feature_sig_unlabel_merged.RData'))
  }


  cat('\n')
  message(crayon::blue('2. Find pairs from these enriched features ... \n'))

  isotope_label_matched <- pbapply::pbmapply(function(x, y){
    cat(x, '@', y, '\n')
    find_pair_feature(target_mz = x,
                      target_rt = y,
                      labeled_feature_table = feature_sig_label,
                      mz_tol = 15,
                      rt_tol = 0.05)
  },
  x = feature_sig_unlabel$mz,
  y = feature_sig_unlabel$rt,
  SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()

  pair_table <- feature_sig_unlabel %>%
    # dplyr::bind_cols(isotope_label_matched) %>%
    dplyr::left_join(isotope_label_matched, by = c('id' = 'source_feature'), multiple = "all") %>%
    dplyr::select(id, matched_feature, mz:rt, actual_mz, actual_rt, mass_shift_label, p_values, q_values, fold_change, avg_abundance_case) %>%
    dplyr::rename('unlabeled_feature_id' = 'id',
                  'labelded_feature_id' = 'matched_feature',
                  'unlabeled_mz' = 'mz',
                  'unlabeled_rt' = 'rt',
                  'labeled_mz' = 'actual_mz',
                  'labeled_rt' = 'actual_rt',
                  'p_value_unlabel' = 'p_values',
                  'q_value_unlabel' = 'q_values',
                  'fold_change_unlabel' = 'fold_change',
                  'average_abundance_unlabel' = 'avg_abundance_case') %>%
    dplyr::filter(!is.na(labelded_feature_id))

  temp_idx <- match(pair_table$labelded_feature_id, raw_data_label$id)
  pair_table <- pair_table %>%
    dplyr::mutate(p_value_label = raw_data_label$p_values[temp_idx],
                  q_value_label = raw_data_label$q_values[temp_idx],
                  fold_change_label = raw_data_label$fold_change[temp_idx],
                  average_abundance_label = raw_data_label$avg_abundance_case[temp_idx])


  message(crayon::yellow('Result: There are ', nrow(pair_table), 'labeled pairs are found.\n'))
  print(knitr::kable(pair_table))

  if (is_recognize_adducts) {
    table_export <- list('raw_data_unlabel' = raw_data_unlabel,
                         'raw_data_label' = raw_data_label,
                         'recognized_peaks_unlabel' = result_peak_recognization_unlabel,
                         'paired_table' = pair_table)
    dir.create(file.path(path, '00_tracer_result'), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(table_export,
                        path = file.path(path, '00_tracer_result', 'tracer_pair_result.xlsx'),
                        format_headers = FALSE)
  } else {
    table_export <- list('raw_data_unlabel' = raw_data_unlabel,
                         'raw_data_label' = raw_data_label,
                         # 'recognized_peaks_unlabel' = result_peak_recognization_unlabel,
                         'paired_table' = pair_table)
    dir.create(file.path(path, '00_tracer_result'), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(table_export,
                        path = file.path(path, '00_tracer_result', 'tracer_pair_result.xlsx'),
                        format_headers = FALSE)
  }


  cat('\n')
  message(crayon::blue('3. Plot these pairs ... \n'))
  message(crayon::blue('3.1 Extract EICs for pairs ...\n'))

  if (case_group[1] %in% list.files(path)) {
    temp_files <- list.files(path = file.path(path, case_group[1]), recursive = TRUE)
    eic_data_unlabel <- extract_eic_data(path = path,
                                         files = file.path(case_group[1], temp_files),
                                         mz_list = pair_table$unlabeled_mz,
                                         mz_tol = mz_tol)
  } else {
    eic_data_unlabel <- extract_eic_data(path = path,
                                         files_pattern = paste0(case_group[1], '.+\\.mzML'),
                                         mz_list = pair_table$unlabeled_mz,
                                         mz_tol = mz_tol)
  }

  if (case_group[2] %in% list.files(path)) {
    temp_files <- list.files(path = file.path(path, case_group[2]), recursive = TRUE)
    eic_data_label <- extract_eic_data(path = path,
                                         files = file.path(case_group[2], temp_files),
                                         mz_list = pair_table$unlabeled_mz,
                                         mz_tol = mz_tol)
  } else {
    eic_data_label <- extract_eic_data(path = path,
                                       files_pattern = paste0(case_group[2], '.+\\.mzML'),
                                       mz_list = pair_table$labeled_mz,
                                       mz_tol = mz_tol)
  }


  cat('\n')
  message(crayon::blue('3.2 Plot pairs\n'))

  temp_plot <- plot_isotope_pairs(signif_feature_unlabeled = feature_sig_unlabel,
                                  signif_feature_labeled = feature_sig_label,
                                  pair_table = pair_table,
                                  eic_data_unlabel = eic_data_unlabel,
                                  eic_data_label = eic_data_label,
                                  mz_tol = mz_tol)

  dir.create(file.path(path, '00_tracer_result'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(plot = temp_plot,
                  filename = file.path(path, '00_tracer_result', 'isotope_pair_plot.pdf'),
                  width = 8, height = 6)

  message(crayon::blue('Done!\n'))


}



################################################################################
# modify_msdial_table ----------------------------------------------------------

# peak_table_unlabel <- '230830_Csp_YgeX_WT_UA.xlsx'
# peak_table_label <- '230830_Csp_YgeX_WT_13CUA.xlsx'
#
# raw_table <- modify_msdial_table(peak_table = '230830_Csp_YgeX_WT_UA.xlsx',
#                                  path = '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/',
#                                  control_group = "WT_UA",
#                                  case_group = 'ygeX_UA')


# peak_table <- 'ssnA_UA_48h_area.txt'
# peak_table_label <- 'ssnA_13CUA_48h_area.txt'
# path = '~/Project/00_Uric_Acid_project/Data/230928_isotope_tracing/ssnA/'
# control_group = "WT_UA"
# case_group = "ssnA_UA"
#
# raw_table <- modify_msdial_table(peak_table = '230830_Csp_YgeX_WT_UA.xlsx',
#                                  path = '~/Project/00_Uric_Acid_project/20230831_isotope_tracing/',
#                                  control_group = "WT_UA",
#                                  case_group = 'ygeX_UA')



modify_msdial_table <- function(peak_table,
                                path = '.',
                                control_group = "WT_UA",
                                case_group = 'ygeX_UA') {
  if (stringr::str_detect(peak_table, '\\.xlsx')) {
    raw_table <- readxl::read_xlsx(file.path(path, peak_table))
  } else {
    raw_table <- readr::read_tsv(file.path(path, peak_table))
  }

  temp_col_name <- raw_table %>% dplyr::slice(4) %>% as.character()
  if (!any(stringr::str_detect(temp_col_name, control_group))) {
    stop('No control group existed, please check it')
  }
  if (!any(stringr::str_detect(temp_col_name, case_group))) {
    stop('No case group existed, please check it')
  }

  raw_table <- raw_table %>%
    dplyr::slice(-c(1:4))
  colnames(raw_table) <- temp_col_name

  raw_table <- raw_table %>%
    dplyr::select(`Average Mz`,
                  `Average Rt(min)`,
                  `MS1 isotopic spectrum`,
                  `S/N average`,
                  dplyr::starts_with(control_group),
                  dplyr::starts_with(case_group)) %>%
    dplyr::rename('mz' = `Average Mz`,
                  'rt' = `Average Rt(min)`,
                  'ms1_isotopes' = `MS1 isotopic spectrum`,
                  'sn' = `S/N average`) %>%
    dplyr::mutate(dplyr::across(-'ms1_isotopes', as.numeric)) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 4)) %>%
    dplyr::mutate(id = paste0(mz, '@', rt)) %>%
    dplyr::select(id, dplyr::everything())

  temp_data <- raw_table %>%
    dplyr::select(dplyr::starts_with(c(control_group, case_group)))
  idx_control <- colnames(temp_data) %>% stringr::str_detect(pattern = control_group) %>% which()
  idx_case <- colnames(temp_data) %>% stringr::str_detect(pattern = case_group) %>% which()
  p_values <- apply(temp_data, 1, function(x){
    result <- t.test(x[idx_control], x[idx_case], alternative = 'two.sided', var.equal = TRUE)
    # result <- wilcox.test(x[idx_control], x[idx_case])
    result$p.value
  })
  q_values <- p.adjust(p_values, method = 'BH')

  fold_change <- apply(temp_data, 1, function(x){
    mean(x[idx_case])/mean(x[idx_control])
  })

  # average_abundance = apply(temp_data, 1, mean)
  average_abundance_control <- apply(temp_data, 1, function(x){mean(x[idx_control])})
  average_abundance_case <- apply(temp_data, 1, function(x){mean(x[idx_case])})

  result_table <- raw_table %>%
    dplyr::mutate(p_values = p_values,
                  q_values = q_values,
                  fold_change = fold_change,
                  avg_abundance_control = average_abundance_control,
                  avg_abundance_case = average_abundance_case)

  return(result_table)
}


# ################################################################################
# modify_msdial_table_timepoints -----------------------------------------------



modify_msdial_table_timepoints <- function(peak_table,
                                           sample_info = 'sample_info.xlsx',
                                           path = '.',
                                           control_group = "WT",
                                           case_group = 'ygeY',
                                           iso_condition = c('12C', '13C')) {
  iso_condition <- match.arg(iso_condition)
  if (stringr::str_detect(peak_table, '\\.xlsx')) {
    raw_table <- readxl::read_xlsx(file.path(path, peak_table))
  } else {
    raw_table <- readr::read_tsv(file.path(path, peak_table))
  }

  temp_col_name <- raw_table %>% dplyr::slice(4) %>% as.character()
  if (!any(stringr::str_detect(temp_col_name, control_group))) {
    stop('No control group existed, please check it')
  }
  if (!any(stringr::str_detect(temp_col_name, case_group))) {
    stop('No case group existed, please check it')
  }

  raw_table <- raw_table %>%
    dplyr::slice(-c(1:4))
  colnames(raw_table) <- temp_col_name


  # check raw table and sample info
  sample_info <- readxl::read_xlsx(file.path(path, sample_info))
  # check groups whether included in the sample info
  if (!all(c(control_group, case_group) %in% sample_info$group)) {
    stop('The control_group and case_group do not included in the sample info file\n')
  }

  # check samples whether included in the sample info
  sample_info <- sample_info %>%
    dplyr::filter(type == 'sample') %>%
    dplyr::filter(stable_isotope_labeling == iso_condition) %>%
    dplyr::arrange(match(group, c(control_group, case_group), timepoint, sample_name))

  temp_sample_name <- sample_info %>%
    dplyr::pull(sample_name)

  if(!all(temp_sample_name %in% colnames(raw_table))) {
    temp <- which(!(temp_sample_name %in% colnames(raw_table))) %>%
      temp_sample_name[.]

    stop(length(temp), 'samples are not included in the sample info table\n',
         paste(temp, collapse = ', '))
  }

  raw_table <- raw_table %>%
    dplyr::select(`Average Mz`,
                  `Average Rt(min)`,
                  `MS1 isotopic spectrum`,
                  `S/N average`,
                  dplyr::all_of(temp_sample_name)) %>%
    dplyr::rename('mz' = `Average Mz`,
                  'rt' = `Average Rt(min)`,
                  'ms1_isotopes' = `MS1 isotopic spectrum`,
                  'sn' = `S/N average`) %>%
    dplyr::mutate(dplyr::across(-'ms1_isotopes', as.numeric)) %>%
    dplyr::mutate(mz = round(mz, 4),
                  rt = round(rt, 4)) %>%
    dplyr::mutate(id = paste0(mz, '@', rt)) %>%
    dplyr::select(id, dplyr::everything())

  # extract peak abundance as matrix
  temp_data <- raw_table %>%
    dplyr::select(dplyr::all_of(temp_sample_name)) %>%
    as.matrix()

  # state for each time
  stat_result <- lapply(unique(sample_info$timepoint), function(tp){
    idx_control <- sample_info %>%
      dplyr::filter(timepoint == tp & group == control_group) %>%
      dplyr::pull(sample_name) %>%
      match(., colnames(temp_data))
    idx_case <- sample_info %>%
      dplyr::filter(timepoint == tp & group == case_group) %>%
      dplyr::pull(sample_name) %>%
      match(., colnames(temp_data))


    p_values <- apply(temp_data, 1, function(x){
      if (sd(x[c(idx_control, idx_case)]) == 0) {
        return(1)
      } else if (sd(x[idx_control]) == 0 & sd(x[idx_case]) == 0) {
        return(0)
      } else {
        result <- t.test(x[idx_control], x[idx_case], alternative = 'two.sided', var.equal = TRUE)
        return(result$p.value)
      }
    })

    q_values <- p.adjust(p_values, method = 'BH')

    fold_change <- apply(temp_data, 1, function(x){
      mean(x[idx_case])/mean(x[idx_control])
    })

    average_abundance_control <- apply(temp_data, 1, function(x){mean(x[idx_control])})
    average_abundance_case <- apply(temp_data, 1, function(x){mean(x[idx_case])})

    result <- list(p_values = p_values,
                   q_values = q_values,
                   fold_change = fold_change,
                   average_abundance_control = average_abundance_control,
                   average_abundance_case = average_abundance_case)
  })

  stat_result <- stat_result %>% bind_cols()
  colnames(stat_result) <- c('p_values', 'q_values', 'fold_change', paste0('avg_', control_group), paste0('avg_', case_group)) %>%
    rep(., times = length(unique(sample_info$timepoint))) %>%
    paste0(., '_T', rep(unique(sample_info$timepoint), each = 5))

  # add stat result to raw_table
  result_table <- raw_table %>%
    dplyr::bind_cols(stat_result)

  # if there is multiple timepoint
  if (length(unique(sample_info$timepoint)) > 1) {
    # test slope difference
    temp_data2 <- temp_data %>%
      tibble::as_tibble() %>%
      dplyr::mutate(id = raw_table$id) %>%
      tidyr::pivot_longer(cols = -id, names_to = 'sample_name', values_to = 'abundance')

    temp_idx <- match(temp_data2$sample_name, sample_info$sample_name)
    temp_data2 <- temp_data2 %>%
      dplyr::mutate(group = sample_info$group[temp_idx],
                    timepoint = sample_info$timepoint[temp_idx])

    p_value_timepoint <- sapply(unique(temp_data2$id), function(x){
      lm_data <- temp_data2 %>%
        dplyr::filter(id == x)
      lm_model <- lm(abundance ~ timepoint * group, data = lm_data)
      interaction_test <- lmtest::coeftest(lm_model, vcov = sandwich::vcovHC(lm_model, type = "HC1"))
      p_value <- interaction_test[4,4]
      return(p_value)
    })

    q_value_timepoint <- p.adjust(p_value_timepoint, method = 'BH')

    result_table <- result_table %>%
      dplyr::mutate(p_value_timepoint = p_value_timepoint,
                    q_value_timepoint = q_value_timepoint)


    # monotonicity check, only check case group
    temp_data3 <- result_table %>%
      dplyr::select(
        # id,
        # dplyr::starts_with(paste0('avg_', control_group)),
        dplyr::starts_with(paste0('avg_', case_group))) %>%
      as.matrix()


    monotonicity <- apply(temp_data3, 1, function(x){
      cor_spearman <- cor(x, unique(sample_info$timepoint), method = "spearman")
    })

    result_table <- result_table %>%
      dplyr::mutate(monotonicity = monotonicity)
  } else {
    result_table <- result_table %>%
      dplyr::mutate(p_value_timepoint = 1, q_value_timepoint = 1, monotonicity = 0)
  }

  return(result_table)
}



################################################################################
# find_pair_feature ------------------------------------------------------------

# target_mz = 105.0645
# target_rt = 7.629
# mz_tol <- 15
# rt_tol <- 0.05
# labeled_feature_table = feature_sig_label

find_pair_feature <- function(target_mz,
                              target_rt,
                              labeled_feature_table,
                              mz_tol,
                              rt_tol) {
  # browser()

  # predict target_mz formula
  pred_formula_table <- MassToolsMjhelf::calcMF(mz = target_mz, z = 1, ppm = mz_tol)

  # if the formula can't be predicted, export the NA
  if (length(pred_formula_table) < 1)  {
    output <- data.frame(source_feature = NA,
                         matched_feature = NA,
                         mass_shift_label = NA,
                         theo_mz = NA,
                         actual_mz = NA,
                         mz_error = NA,
                         theo_rt = NA,
                         actual_rt = NA,
                         rt_error = NA)
    return(output)
  }

  carbon_number_max <- pred_formula_table$MF %>%
    lapply(function(x){
      temp_result <- rcdk::get.formula(x)
      temp_result@isotopes %>%
        as_tibble() %>%
        dplyr::filter(isoto == 'C') %>%
        dplyr::pull(number) %>%
        as.numeric()
    }) %>%
    unlist() %>%
    unique() %>%
    max()

  possible_mz <- target_mz + c(1:carbon_number_max)*1.00335

  possible_table <- data.frame(mz = possible_mz,
                               rt = target_rt,
                               mass_shift_label = paste0('13C*', c(1:carbon_number_max)),
                               stringsAsFactors = FALSE)

  temp_label_table <- labeled_feature_table %>%
    dplyr::select(mz, rt) %>%
    as.data.frame()

  matched_result <- masstools::mz_rt_match(
    possible_table,
    temp_label_table,
    mz.tol = mz_tol,
    rt.tol = rt_tol,
    rt.error.type = 'abs'
  )

  if (is.null(matched_result)) {
    output <- data.frame(source_feature = NA,
                         matched_feature = NA,
                         mass_shift_label = NA,
                         theo_mz = NA,
                         actual_mz = NA,
                         mz_error = NA,
                         theo_rt = NA,
                         actual_rt = NA,
                         rt_error = NA)
    return(output)
  }

  # if (nrow(matched_result) == 0) {
  #   return(NA)
  # }

  # index
  idx1 <- matched_result$Index1
  idx2 <- matched_result$Index2

  # export matched pair table
  matched_feature <- labeled_feature_table[idx2,] %>%
    dplyr::pull(id)

  source_feature <- paste0(target_mz, '@', target_rt)

  label_isotope <- possible_table[idx1,] %>% dplyr::pull(mass_shift_label)

  output <- matched_result %>%
    dplyr::mutate(source_feature = source_feature,
                  matched_feature = matched_feature,
                  mass_shift_label = label_isotope) %>%
    dplyr::select(source_feature, matched_feature, mass_shift_label, mz1, mz2, `mz error`, rt1, rt2, `rt error`) %>%
    dplyr::rename('theo_mz' = 'mz1',
                  'actual_mz' = 'mz2',
                  'mz_error' = 'mz error',
                  'theo_rt' = 'rt1',
                  'actual_rt' = 'rt2',
                  'rt_error' = 'rt error')

  return(output)
}





################################################################################
# extract_eic_data -------------------------------------------------------------


# path <- '.'
# files_pattern <- 'ygeX_UA.+\\.mzML'
# mz_list <- pair_table$unlabeled_mz
# mz_tol = 10

extract_eic_data <- function(path,
                             files = NULL,
                             files_pattern = NULL, # 'ygeX_UA.+\\.mzML'
                             mz_list,
                             mz_tol = 10) {
  # browser()
  require(data.table)
  require(RaMS)

  if (is.null(files)) {
    if (is.null(files_pattern)) {
      stop('Please provide files or files_pattern\n')
    }

    files <- list.files(path = path, pattern = files_pattern, recursive = TRUE)

    if (length(files) == 0) {
      stop('No effective mzML data found\n')
    }
  }

  msdata <- RaMS::grabMSdata(file.path(path, files), grab_what = 'MS1')

  # unlabeled EIC
  mz_list <- unique(round(mz_list, 4))

  temp_mz_tol <- lapply(mz_list, function(x){
    RaMS::pmppm(x, ppm = mz_tol)
  })

  idx_list <- lapply(temp_mz_tol, function(x){
    which(msdata$MS1$mz >= x[1] & msdata$MS1$mz <= x[2])
  })

  eic_data_export <- mapply(function(x, y){
    msdata$MS1[x,] %>%
      dplyr::mutate(label = stringr::str_replace(filename, '\\_\\d+\\.mzML', ''),
                    precursor_mz = as.character(y))
  },
  x = idx_list,
  y = mz_list,
  SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()

  return(eic_data_export)
}

################################################################################
# plot_isotope_pairs -----------------------------------------------------------

plot_isotope_pairs <- function(signif_feature_unlabeled,
                               signif_feature_labeled,
                               pair_table,
                               eic_data_unlabel,
                               eic_data_label,
                               col_abundance = 'avg_abundance_case',
                               mz_tol = 10) {

  # point table
  # temp_labeled <- signif_feature_unlabeled %>%
  #   dplyr::select(id, mz, rt, get(average_abundance), q_value, fold_change) %>%
  #   dplyr::mutate(group_label = 'Unlabeled',
  #                 delta_mz = 0.1)
  #
  # temp_unlabeled <- signif_feature_labeled %>%
  #   dplyr::select(id, mz, rt, average_abundance, q_value, fold_change) %>%
  #   dplyr::mutate(group_label = 'Labeled',
  #                 delta_mz = -0.1)

  temp_labeled <- signif_feature_unlabeled %>%
    dplyr::select(id, mz, rt, col_abundance, q_values, fold_change) %>%
    # dplyr::rename('average_abundance' = col_abundance) %>%
    dplyr::mutate(group_label = 'Unlabeled',
                  delta_mz = 0.1)

  temp_unlabeled <- signif_feature_labeled %>%
    dplyr::select(id, mz, rt, col_abundance, q_values, fold_change) %>%
    dplyr::mutate(group_label = 'Labeled',
                  delta_mz = -0.1)

  point_table <- temp_unlabeled %>%
    dplyr::bind_rows(temp_labeled)

  # pair table
  pair_table <- pair_table %>%
    dplyr::mutate(avg_mz = (unlabeled_mz + labeled_mz)/2) %>%
    dplyr::mutate(mz_adj1 = unlabeled_mz - avg_mz,
                  mz_adj2 = labeled_mz - avg_mz)

  # replace delta_mz in point table
  temp_idx1 <- match(pair_table$unlabeled_feature_id, point_table$id)
  temp_idx2 <- match(pair_table$labelded_feature_id, point_table$id)

  point_table$delta_mz[temp_idx1] <- pair_table$mz_adj1
  point_table$delta_mz[temp_idx2] <- pair_table$mz_adj2

  # segment table1
  segment_table <- pair_table %>%
    dplyr::mutate(rt_mean = (unlabeled_rt + labeled_rt)/2) %>%
    dplyr::mutate(rt_adj = rt_mean - 0.1) %>%
    dplyr::select(mz_adj1, mz_adj2, rt_mean, rt_adj, mass_shift_label)

  # rescale EICs
  temp_scale <- round(pair_table$labeled_mz - pair_table$unlabeled_mz) %>% max()
  eic_data_unlabel <- eic_data_unlabel %>%
    dplyr::mutate(int_adj = int/max(int)*(-temp_scale/2) - 0.1)
  eic_data_label <- eic_data_label %>%
    dplyr::mutate(int_adj = int/max(int)*(temp_scale/2) + 0.1)

  temp_plot <- ggplot2::ggplot(point_table) +
    ggplot2:: geom_line(data = eic_data_label,
                        ggplot2::aes(x = rt, y = int_adj, color = precursor_mz),
                        alpha = 0.6) +
    ggplot2::geom_line(data = eic_data_unlabel,
                       ggplot2::aes(x = rt, y = int_adj, color = precursor_mz), alpha = 0.6) +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj1,
                                       xend = rt_adj,
                                       yend = mz_adj2),
                          color = 'gray') +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj1,
                                       xend = rt_mean,
                                       yend = mz_adj1),
                          color = 'gray') +
    ggplot2::geom_segment(data = segment_table,
                          ggplot2::aes(x = rt_adj,
                                       y = mz_adj2,
                                       xend = rt_mean,
                                       yend = mz_adj2),
                          color = 'gray') +
    ggplot2::geom_point(ggplot2::aes(x = rt,
                                     y = delta_mz,
                                     fill = group_label,
                                     size = fold_change),
                        alpha = 0.5, shape = 21) +
    ggrepel::geom_text_repel(
      data = subset(point_table, abs(delta_mz) > 0.5 ),
      ggplot2::aes(x = rt, y = delta_mz, label = id),
      position = position_jitter(seed = 1)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_fill_manual(values = c('Labeled' = 'red',
                                          'Unlabeled' = 'black'),
                               label = c('Labeled' = 'C13-UA',
                                         'Unlabeled' = 'UA'),
                               name = 'Group') +
    ggplot2::scale_colour_discrete(name = 'Precursor m/z') +
    ggplot2::scale_size_continuous(range = c(1, 10),
                                   name = paste0('Fold_change (Mutant/WT)')) +
    ggplot2::scale_y_continuous(breaks = c((seq(-floor(temp_scale/2), -1, by = 1)) - 0.1,
                                           0,
                                           (seq(1, floor(temp_scale/2), by = 1)) + 0.1),
                                labels = seq(-floor(temp_scale/2), floor(temp_scale/2), by = 1)) +
    ggplot2::ylab('m/z - average m/z') +
    ggplot2::xlab('RT (min)') +
    ZZWTheme() +
    ggplot2::theme(legend.position = 'right')

  return(temp_plot)

}


################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
Version 0.1.6
-------------
Authors: Zhiwei Zhou
Maintainer: Zhiwei Zhou

Updates
-------------
o Fixed bug: Assign ppc score as 0 if no EIC of base peak can be extracted
o Add base peak average intensity information in the peak group, then modify the order of merge


v0.1.3
-------------
o Fixed bug: Add a function for multiple time point data analysis

v0.1.4
-------------
o Set an option to skip ppc calculation

v0.1.5
-------------
o Correct 'Average_Abundance' name

v0.1.6
-------------
o Fix bug: self_check_isf
o Adjust ppc calculation (5 points needed)
o Adjust the formula prediction threshold (ppm: 10 -> 15)

v0.1.7
-------------
o Add README.md
o Commit and submit to GitHub
")
}
